#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Source: https://github.com/klamkiew/viralclust/commit/3203e6de334c6834877dbdffecff70df07ed80d7
# "Clustering of viral genomes based on different algorithms and metrices"
# (Author: kevin.lamkiewicz@uni-jena.de)

"""
Usage:
  hdbscan_virus.py [options] <inputSequences> <lineageDict>

Options:
  -h, --help                              Show this help message and exits.
  -v, --verbose                           Get some extra information from viralClust during calculation. [Default: False]
  -o DIR, --output DIR                    Specifies the output directory of viralClust. [Default: pwd]
  -p PROCESSES, --process PROCESSES       Specify the number of CPU cores that are used. [Default: 1]
  -k KMER, --kmer KMER                    Length of the considered kmer. [Default: 5]
  --umap_metric UMAP_METRIC                         [Default: cosine] 
  --n_components N_COMPONENTS                         [Default: 2] 
  --n_neighbors N_NEIGHBORS                          [Default: 15]
  --min_dist MIN_DIST                              [Default: 0.1]
  --hdbscan_metric HDBSCAN_METRIC                       [Default: cosine] 
  --min_samples MIN_SAMPLES                          [Default: 5] 
  --min_cluster_size MIN_CLUSTER_SIZE                      [Default: 5]  
  --cluster_selection_epsilon CLUSTER_SELECTION_EPSILON             [Default: 0.05] 

"""

####################################################################
# IMPORTS 
####################################################################

import sys
import os
import logging
import glob
import shutil

import numpy as np
import pandas as pd
from multiprocessing import Pool

from docopt import docopt

from collections import Counter
import multiprocessing as mp
import itertools
import math
import subprocess

import scipy
import random
import umap.umap_ as umap
import umap.plot
import matplotlib
matplotlib.use("Agg") # use non-interactive backend to not generate figure windows when already being disconnected from the computing server
import matplotlib.pyplot as plt
import hdbscan
from sklearn.preprocessing import normalize
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score


########################################################
# DATA STRUCTURES
########################################################
class Clusterer(object):
  """
  """

  id2header = {}
  d_profiles = {}
  header2id = {}
  header2lineage = {}
  #header2primer = {}
  dim = 0
  matrix = np.empty(shape=(dim,dim))

  genomeOfInterest = ''

  scipyDistances = {
      'euclidean' : scipy.spatial.distance.euclidean ,
      'manhatten' : scipy.spatial.distance.cityblock ,
      'chebyshev' : scipy.spatial.distance.chebyshev  ,
      'minkwoski': scipy.spatial.distance.minkowski ,
      'canberra' : scipy.spatial.distance.canberra ,
      'braycurtis' : scipy.spatial.distance.braycurtis ,
      'mahalanobis' : scipy.spatial.distance.mahalanobis ,
      'wminkowski' : scipy.spatial.distance.wminkowski ,
      'seuclidean' : scipy.spatial.distance.seuclidean ,
      'cosine' : scipy.spatial.distance.cosine
  }

  def __init__(self, logger, sequenceFile, lineageDict, outdir,  k, proc, umap_metric, n_components, n_neighbors,  min_dist, hdbscan_metric, min_cluster_size, min_samples, cluster_selection_epsilon):
    """
    """

    self.lineageDict = lineageDict
    #self.primerDict = primerDict
    self.sequenceFile = sequenceFile
    self.outdir = outdir

    self.reducedSequences = f"{self.outdir}/{os.path.basename(sequenceFile)}"

    self.k = k
    nucleotides = set(["A","C","G","T"])
    self.allKmers = {''.join(kmer):x for x,kmer in enumerate(itertools.product(nucleotides, repeat=self.k))}

    self.proc = proc
    self.umap_metric = umap_metric
    self.hdbscan_metric = hdbscan_metric
    self.n_neighbors = n_neighbors
    self.min_dist = min_dist
    self.n_components = n_components
    self.min_cluster_size = min_cluster_size
    self.min_samples = min_samples
    self.cluster_selection_epsilon = cluster_selection_epsilon

    self.d_sequences = {}
    self.allCluster = []
    self.clusterlabel = []
    self.probabilities = []


  def __parse_fasta(self, filePath):
    """
    """
    with open(filePath, 'r') as inputStream:
      header = ''
      seq = ''
      for line in inputStream:
        if line.startswith(">"):
          if header:
            yield (header, seq)
          header = line.rstrip("\n").replace(':','_').replace(' ','_').strip(">_")
          seq = ''
        else:
          seq += line.rstrip("\n").upper().replace('U','T')
      yield (header, seq)

  def read_sequences(self):
    """
    """
    idHead = -1
    fastaContent = {}
    for header, sequence in self.__parse_fasta(self.sequenceFile):
      idHead += 1
      if header in Clusterer.header2id.keys():
        logger.warn(f"{header} already stored with id {Clusterer.header2id[header]}")
      Clusterer.id2header[idHead] = header
      Clusterer.header2id[header] = idHead
      fastaContent[idHead] = sequence
    # print(f"size id2header and header2id: {len(Clusterer.id2header), len(Clusterer.header2id)}")
    # print(f"unique values id2header and header2id: {len(set(list(Clusterer.id2header.values()))), len(set(list(Clusterer.header2id.values())))}")

    Clusterer.dim = len(fastaContent)
    Clusterer.matrix = np.zeros(shape=(Clusterer.dim, Clusterer.dim), dtype=float)
    # for hdbscan, the data set size has to be larger or equal to min_samples
    if len(fastaContent) < 21:
      return 1
    return fastaContent


  def profile(self, entry):
    """
    """
    header, sequence = entry
    profile = [0]*len(self.allKmers)
    for k in iter([sequence[start : start + self.k] for start in range(len(sequence) - self.k)]):
        try:
          profile[self.allKmers[k]] += 1
        except KeyError:
          continue
    kmerSum = sum(profile)
    try:
      profile = list(map(lambda x: x/kmerSum, profile))
      return (header, profile)
    except ZeroDivisionError:
      logger.warn(Clusterer.id2header[header] + " skipped, due to too many N's.")
      return(None, None)


  def determine_profile(self, proc):

    self.d_sequences = self.read_sequences()
    if self.d_sequences == 1:
      logger.warn(f"{self.d_sequences} sequence(s) are too few for clustering")
      import shutil
      shutil.copyfile(self.sequenceFile, f'{self.outdir}/{os.path.splitext(os.path.basename(self.sequenceFile))[0]}_hdbscan.fasta')
      return 1

    p = Pool(self.proc)
    allProfiles = p.map(self.profile, self.d_sequences.items())
    p.close()
    p.join()
    for header, profile in allProfiles:
      #if header:
      Clusterer.d_profiles[header] = profile

  def get_lineage_annotations(self):
    df = pd.read_csv(self.lineageDict, sep=',')
    data = {r.header: r.lineage for i,r in df.iterrows()}
    for index, profile in Clusterer.d_profiles.items():
      header = Clusterer.id2header[index]
      if 'sub' in header:
        Clusterer.header2lineage[header] = data['_'.join(header.split('_')[:-1])]
      else:  
        Clusterer.header2lineage[header] = data[header]
      # Clusterer.header2primer[header] = primer_data[header]
    # print("Number of profiles, number of lineage assignments and primer locations:")
    # print(len(Clusterer.d_profiles), len(Clusterer.header2lineage), len(Clusterer.header2primer))
    # print(f"Number of profiles, number of lineage assignments: {len(Clusterer.d_profiles)}, {len(Clusterer.header2lineage)}")



  def apply_umap(self):
    """
    """
    profiles = [(idx,profile) for idx, profile in Clusterer.d_profiles.items() if idx in self.d_sequences]
    vector = np.array([x[1] for x in profiles])
    clusterable_embedding = umap.UMAP(
                n_neighbors=2*self.n_neighbors,
                min_dist=0,
                n_components=self.n_components,
                random_state=42,
                metric=self.umap_metric,
                unique=True
            ).fit_transform(vector)
    visual_embedding = umap.UMAP(
                n_neighbors=self.n_neighbors,
                min_dist=self.min_dist,
                n_components=self.n_components,
                random_state=42,
                metric=self.umap_metric,
                unique=True
            ).fit(vector)

    # Plot UMAPped reads with primer location
    # fig, ax = plt.subplots()
    # ax = umap.plot.points(visual_embedding, labels=np.array(list(Clusterer.header2primer.values())))
    # ax.set_title('UMAP embedding of input reads (colored by primer location)')
    # plt.savefig('umap_primerplot.png')
    # plt.close()

    clusterer = hdbscan.HDBSCAN(min_cluster_size=self.min_cluster_size, min_samples=self.min_samples, cluster_selection_epsilon=self.cluster_selection_epsilon)
    if self.hdbscan_metric == "cosine":
      clusterable_embedding = normalize(clusterable_embedding,norm='l2')
      clusterer.fit(clusterable_embedding)
    else:
      clusterer.fit(clusterable_embedding,metric=self.hdbscan_metric)

    self.clusterlabel = clusterer.labels_
    self.probabilities = clusterer.probabilities_
    #   if len(set(self.clusterlabel)) == 1:
    #     raise TypeError

    # except TypeError:
    #   import shutil
    #   shutil.copyfile(self.sequenceFile, f'{self.outdir}/{os.path.splitext(os.path.basename(self.sequenceFile))[0]}_hdbscan.fasta')
    #   return 1

    self.allCluster = list(zip([x[0] for x in profiles], self.clusterlabel))


    with open(f'{self.outdir}/cluster.txt', 'w') as outStream:
      for i in set(self.clusterlabel):
        cluster_profiles = []
        with open(f'{self.outdir}/cluster{i}.fasta', 'w') as fastaOut:
          outStream.write(f">Cluster {i}\n")
          for idx, label in self.allCluster:
            if label == i:
              cluster_profiles.extend([t[1] for t in profiles if t[0]==idx])
              outStream.write(f"{Clusterer.id2header[idx]}\n")
              fastaOut.write(f">{Clusterer.id2header[idx]}\n{self.d_sequences[idx].split('X'*10)[0]}\n")
        outStream.write("\n")
        print(f"Cluster {i}:\t{np.array(cluster_profiles).shape[0]} reads")
        cluster_profile_mean = np.mean(np.array(cluster_profiles), axis=0)
        plt.subplots(figsize=(12,5))
        plt.hist(cluster_profile_mean)
        plt.title(f'Mean Kmer Histogram Cluster {i}', fontsize=13)
        plt.savefig(f'mean-histogram_cluster_{i}.png')
        #plt.close()

    # # Plot clustering result in UMAP space
    # fig, ax =plt.subplots()
    # ax.scatter(clusterable_embedding[:,0], clusterable_embedding[:,1], c=self.clusterlabel, label=self.clusterlabel, cmap='Paired')
    # ax.legend()
    # ax.set_title('Normalized UMAP embedding of input reads (colored by HDBSCAN cluster)')
    # plt.savefig('normalized_umap_clusterplot.png')
    # plt.close()
    plt.subplots()
    umap.plot.points(visual_embedding, labels=self.clusterlabel)
    plt.title(f'UMAP embedding colored by HDBSCAN cluster')
    plt.savefig('umap_clusterplot.png')
    #plt.close()
    if 'dummy' not in self.lineageDict:
      #print(f"Size of UMAP clusterable data, visual data, and number of lineage labels:, {len(clusterable_embedding)}, {len(visual_embedding.embedding_)}, {len(Clusterer.header2lineage.values())}")
      ars = adjusted_rand_score(np.array(list(Clusterer.header2lineage.values())), self.clusterlabel)
      mi = adjusted_mutual_info_score(np.array(list(Clusterer.header2lineage.values())), self.clusterlabel)
      print(f"Adjusted rand score: {ars}\nAdjusted mutual information: {mi}")
      # Plot UMAPped reads with lineage assignment
      plt.subplots()
      umap.plot.points(visual_embedding, labels=np.array(list(Clusterer.header2lineage.values())))
      plt.title('UMAP embedding of input reads (colored by lineage assignment)')
      plt.savefig('umap_lineageplot.png')
      #plt.close()



########################################################
# CLUSTER PIPELINE
########################################################
def __abort_cluster(clusterObject, filename):
    logger.warn(f"Too few sequences for clustering in {os.path.basename(filename)}. No subcluster will be created.")
    del clusterObject

def perform_clustering():

  multiPool = Pool(processes=proc)
  virusClusterer = Clusterer(logger, inputSequences, lineageDict, outdir,  k, proc, umap_metric, n_components, n_neighbors, min_dist, hdbscan_metric, min_cluster_size, min_samples, cluster_selection_epsilon)

  logger.info("Determining k-mer profiles for all sequences.")
  code = virusClusterer.determine_profile(multiPool)
  if code == 1:
    __abort_cluster(virusClusterer, inputSequences)
    return 0
  if 'dummy' not in virusClusterer.lineageDict:
    logger.info("Read lineage annotations for all sequences")
    virusClusterer.get_lineage_annotations()

  logger.info("Clustering with UMAP and HDBSCAN.")
  code = virusClusterer.apply_umap()
  clusterInfo = virusClusterer.clusterlabel
  noise = np.count_nonzero(clusterInfo == -1)
  N_reads = virusClusterer.dim
  logger.info(f"Summarized {N_reads} sequences into {clusterInfo.max()+1} clusters. Filtered {noise} sequences due to uncertainty.")
  logger.info(f"Amount of clustered data: {round((N_reads-noise)/float(N_reads)*100,2)}%")

  
  return 0



#########################################################
# READ INPUT + LOGGER FUNCTIONALITIES
########################################################

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn


def create_logger():
    """
    doc string.
    """

    logger = logging.getLogger()
    handle = logging.StreamHandler()
    logger.addHandler(handle)
    return logger

def create_outdir(outdir):
    try:
      os.makedirs(outdir)
      logger.info(f"Creating output directory: {outdir}")
    except FileExistsError:
      logger.warning(f"!!! The output directory exists. Files will be overwritten !!!.")

def parse_arguments(d_args):
  """
  Parse all given arguments and check for error (e.g. file names).

  Arguments:
  d_args -- dict with input parameters and their values

  Returns:
  parsed and checked parameter list
  """

  verbose = d_args['--verbose']
  if verbose:
    logger.setLevel(logging.INFO)

  inputSequences = d_args['<inputSequences>']
  print('inputSequences', inputSequences)
  if not os.path.isfile(inputSequences):
    logger.error("!!! Couldn't find input sequences. Check your file !!!")
    exit(1)

  lineageDict = d_args['<lineageDict>']
  print('lineageDict', lineageDict)
  if not os.path.isfile(lineageDict):
      logger.error("!!! Couldn't find lineage dictionary. Check your file !!!")
      exit(1)

  # primerDict = d_args['<primerDict>']
  # print('primerDict', primerDict)
  # if not os.path.isfile(primerDict):
  #     logger.error("!!! ouldn't find primer dictionary. Check your file !!!")
  #     exit(1)

  try:
    k = int(d_args['--kmer'])
    print('k', k)
  except ValueError:
    logger.error("!!! Invalid parameter for k-mer size. Please input a number. !!!")
    exit(2)
    
  try:
    proc = int(d_args['--process'])
    print('proc', proc)
  except ValueError:
    logger.error("!!! Invalid number for CPU cores. Please input a number. !!!")
    exit(2)

  output = d_args['--output']
  print('output', output)
  if output == 'pwd':
    output = os.getcwd()
  create_outdir(output)

  METRICES = [
              'euclidean', 'manhatten', 'chebyshev',
              'minkwoski', 'canberra', 'braycurtis',
              'mahalanobis', 'wminkowski',
              'seuclidean', 'cosine'
  ]
  hdbscan_metric = d_args['--hdbscan_metric']
  print('hdbscan_metric', hdbscan_metric)
  if hdbscan_metric not in METRICES:
    log.warn(f"!!! Invalid hdbscan_metric chosen. Will default back to cosine distance. !!!")
    hdbscan_metric = 'cosine'
  
  umap_metric = d_args['--umap_metric']
  print('umap_metric', umap_metric)
  if umap_metric not in METRICES:
    log.warn(f"!!! Invalid umap_metric chosen. Will default back to cosine distance. !!!")
    umap_metric = 'cosine'

  try:
    n_neighbors = int(d_args['--n_neighbors'])
    print('n_neighbors', n_neighbors)
    if n_neighbors <= 0:
      raise ValueError
  except ValueError:
      log.error("!!! Invalid parameter for --n_neighbors. Please input a positive integer. !!!")
      exit(2)

  try:
    min_dist = float(d_args['--min_dist'])
    print('min_dist', min_dist)
    if not 0.0 <= min_dist < 1:
      raise ValueError
  except ValueError:
    log.error("!!! Invalid parameter for --min_dist. Please input a number between [0,1). !!!")
    exit(2)

  try:
    n_components = int(d_args['--n_components'])
    print('n_components', n_components)
    if n_components < 1:
      raise ValueError
  except ValueError:
      log.error("!!! Invalid parameter for --n_components. Please input a positive integer. !!!")
      exit(2)
  
  try:
    min_cluster_size = int(d_args['--min_cluster_size'])
    print('min_cluster_size', min_cluster_size)
    if min_cluster_size < 1:
      raise ValueError
  except ValueError:
      log.error("!!! Invalid parameter for --min_cluster_size. Please input a positive integer. !!!")
      exit(2)

  if d_args['--min_samples'] == "min_cluster_size":
    min_samples = min_cluster_size
    print('min_samples', min_samples)
  else:
    try:
      min_samples = int(d_args['--min_samples'])
      print('min_samples', min_samples)
      if min_samples < 1:
        raise ValueError
    except ValueError:
        log.error("!!! Invalid parameter for --min_samples. Please input a positive integer. !!!")
        exit(2)

  try:
    cluster_selection_epsilon = float(d_args['--cluster_selection_epsilon'])
    print('cluster_selection_epsilon', cluster_selection_epsilon)
    if cluster_selection_epsilon < 0:
      raise ValueError
  except ValueError:
      log.error("!!! Invalid parameter for --cluster_selection_epsilon. Please input a positive integer. !!!")
      exit(2)


  return (inputSequences, lineageDict, output, k, proc, umap_metric, n_components, n_neighbors, min_dist, hdbscan_metric, min_samples, min_cluster_size, cluster_selection_epsilon)



###########################################################################
# MAIN
###########################################################################
if __name__ == "__main__":
  logger = create_logger()
  logger.info('Reading input....')
  (inputSequences, lineageDict,  outdir, k, proc, umap_metric, n_components, n_neighbors, min_dist, hdbscan_metric, min_samples, min_cluster_size, cluster_selection_epsilon) = parse_arguments(docopt(__doc__))

  logger.info("Starting to cluster data....")
  perform_clustering()
  if os.path.islink(f"{os.path.dirname(outdir)}/latest"):
    os.remove(f"{os.path.dirname(outdir)}/latest")
  os.system(f"ln -s {outdir} {os.path.dirname(outdir)}/latest")

