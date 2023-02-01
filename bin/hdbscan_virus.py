#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
ViralClust is a python program that takes several genome sequences
from different viruses as an input.
These will be clustered these sequences into groups (clades) based
on their sequence similarity. For each clade, the centroid sequence is
determined as representative genome, i.e. the sequence with the lowest
distance to all other sequences of this clade.

Python Dependencies:
  docopt
  BioPython
  colorlog
  numpy
  scipy
  umap-learn
  hdbscan

Other Dependencies:
  cd-hit

Contact:
  kevin.lamkiewicz@uni-jena.de


Usage:
  hdbscan_virus.py [options] <inputSequences> <lineageDict> <primerDict>

Options:
  -h, --help                              Show this help message and exits.
  -v, --verbose                           Get some extra information from viralClust during calculation. [Default: False]
  --version                               Prints the version of viralClust and exits.

  -o DIR, --output DIR                    Specifies the output directory of viralClust. [Default: pwd]
  -p PROCESSES, --process PROCESSES       Specify the number of CPU cores that are used. [Default: 1]


  -k KMER, --kmer KMER                    Length of the considered kmer. [Default: 7]

  --metric METRIC                         Distance metric applied by UMAP (if applied) and HDBSCAN.
                                          The following are supported:
                                          'euclidean', 'manhatten', 'chebyshev', 'minkwoski',
                                          'canberra', 'braycurtis',
                                          'mahalanobis', 'wminkowski', 'seuclidean',
                                          'cosine'.
                                          If an invalid metric is set, ViralClust will default back to 
                                          the cosine distance.
                                          [Default: cosine]
  --neighbors NEIGHBORS                   Number of neighbors considered by UMAP to reduce the dimension space.
                                          Low numbers here mean focus on local structures within the data, whereas 
                                          larger numbers may loose fine details. [default: 15]
  --dThreshold dTHRESHOLD                 Sets the threshold for the minimum distance of two points in the low-dimensional space.
                                          Smaller thresholds are recommended for clustering and identifying finer topological structures
                                          in the data. [Default: 0.1]
  --dimension DIMENSION                   UMAP tries to find an embedding for the input data that can be represented by a low-dimensional space.
                                          This parameter tells UMAP how many dimensions should be used for the embedding. Lower numbers may result 
                                          in loss of information, whereas larger numbers will increase the runtime. [Default: 10]

  --clusterSize CLUSTERSIZE               This parameter forces HDBSCAN to form cluster with a size larger-equal to CLUSTERSIZE.
                                          Be aware that some data points (i.e. genomes) or even whole subcluster are considered as noise, if this parameter is set too high.
                                          E.g., if a very distinct viral genus has 40 genomes and the parameter is set to anything >40, HDBSCAN will not form
                                          the genus specific cluster. [Default: 5]
  --minSample MINSAMPLE                   Intuitively, this parameter declares how conservative clustering is performed. Higher values will lead 
                                          to more points considered noise, whereas a low value causes "edge-cases" to be grouped into a cluster.
                                          The default parameter is the same as CLUSTERSIZE. [Default: CLUSTERSIZE]
  --clusterThreshold CLUSTERTHRESHOLD     A distance threshold. Clusters below this value will be merged. [Default: 0.0]
"""

####################################################################
# IMPORTS 
####################################################################

from re import I
import sys
import os
import logging
import glob
import shutil
from tkinter import N

import numpy as np
import json
import pandas as pd
from multiprocessing import Pool

from datetime import datetime

from colorlog import ColoredFormatter
from docopt import docopt
from Bio import Phylo

from collections import Counter
import multiprocessing as mp
import itertools
import math
import subprocess

import scipy
import random

import umap.umap_ as umap
import umap.plot
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

  def __init__(self, logger, sequenceFile, lineageDict, primerDict, outdir,  k, proc, metric, neighbors,  threshold, dimension, clusterSize, minSample):
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
    self.metric = metric
    self.neighbors = neighbors
    self.threshold = threshold
    self.dimension = dimension
    self.clusterSize = clusterSize
    self.minSample = minSample
    self.clusterThreshold = clusterThreshold

    self.d_sequences = {}
    self.centroids = []
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
        logger.warn(f"{header} already stored with id {Clusterer.qheader2id[header]}")
      Clusterer.id2header[idHead] = header
      Clusterer.header2id[header] = idHead
      fastaContent[idHead] = sequence
    print(f"size id2header and header2id: {len(Clusterer.id2header), len(Clusterer.header2id)}")
    print(f"unique values id2header and header2id: {len(set(list(Clusterer.id2header.values()))), len(set(list(Clusterer.header2id.values())))}")

    Clusterer.dim = len(fastaContent)
    Clusterer.matrix = np.zeros(shape=(Clusterer.dim, Clusterer.dim), dtype=float)
    # for hdbscan, the data set size has to be larger or equal to min_samples
    if len(fastaContent) < 5:
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
      logger.warn("{self.d_sequences} are too few for clustering")
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
    # with open(self.primerDict) as json_file:
    #   primer_data = json.load(json_file)
    #print(f"size id2header and header2id: {len(Clusterer.id2header), len(Clusterer.header2id)}")
    for index, profile in Clusterer.d_profiles.items():
      header = Clusterer.id2header[index]
      if 'sub' in header:
        Clusterer.header2lineage[header] = data['_'.join(header.split('_')[:-1])]
      else:  
        Clusterer.header2lineage[header] = data[header]
      # Clusterer.header2primer[header] = primer_data[header]
    # print("Number of profiles, number of lineage assignments and primer locations:")
    # print(len(Clusterer.d_profiles), len(Clusterer.header2lineage), len(Clusterer.header2primer))
    print(f"Number of profiles, number of lineage assignments: {len(Clusterer.d_profiles)}, {len(Clusterer.header2lineage)}")


  def calc_pd(self, seqs):
    """
    """
    for element in seqs:
      try:
        stuff = (element[0], element[1])
      except TypeError:
        return None
    seq1, profile1 = seqs[0]
    seq2, profile2 = seqs[1]
    dFunction = Clusterer.scipyDistances[self.metric]
    distance = dFunction(profile1, profile2)
    return (seq1, seq2, distance)


  def apply_umap(self):
    """
    """
    profiles = [(idx,profile) for idx, profile in Clusterer.d_profiles.items() if idx in self.d_sequences]
    vector = [x[1] for x in profiles]
    clusterable_embedding = umap.UMAP(
                n_neighbors=2*self.neighbors,
                min_dist=self.threshold,
                n_components=self.dimension,
                random_state=42,
                metric=self.metric
            ).fit_transform(vector)
    visual_embedding = umap.UMAP(
                n_neighbors=self.neighbors,
                min_dist=self.threshold,
                n_components=self.dimension,
                random_state=42,
                metric=self.metric
            ).fit(vector)
    print(f"Size of UMAP clusterable data, visual data, and number of lineage labels:, {len(clusterable_embedding)}, {len(visual_embedding.embedding_)}, {len(Clusterer.header2lineage.values())}")
    # Plot UMAPped reads with lineage assignment
    fig, ax = plt.subplots()
    ax = umap.plot.points(visual_embedding, labels=np.array(list(Clusterer.header2lineage.values())))
    ax.set_title('UMAP embedding of input reads (colored by lineage assignment)')
    plt.savefig('umap_lineageplot.png')
    plt.close()

    # Plot UMAPped reads with primer location
    # fig, ax = plt.subplots()
    # ax = umap.plot.points(visual_embedding, labels=np.array(list(Clusterer.header2primer.values())))
    # ax.set_title('UMAP embedding of input reads (colored by primer location)')
    # plt.savefig('umap_primerplot.png')
    # plt.close()

    clusterer = hdbscan.HDBSCAN(min_cluster_size=self.clusterSize, min_samples=self.minSample, cluster_selection_epsilon=self.clusterThreshold)
    if self.metric == "cosine":
      clusterable_embedding = normalize(clusterable_embedding,norm='l2')
      clusterer.fit(clusterable_embedding)
    else:
      clusterer.fit(clusterable_embedding,metric=self.metric)

    self.clusterlabel = clusterer.labels_
    self.probabilities = clusterer.probabilities_
    #   if len(set(self.clusterlabel)) == 1:
    #     raise TypeError

    # except TypeError:
    #   import shutil
    #   shutil.copyfile(self.sequenceFile, f'{self.outdir}/{os.path.splitext(os.path.basename(self.sequenceFile))[0]}_hdbscan.fasta')
    #   return 1

    adj_rs = adjusted_rand_score(list(Clusterer.header2lineage.values()), self.clusterlabel)
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
        fig, ax = plt.subplots(figsize=(12,5))
        ax.hist(cluster_profile_mean)
        plt.title(f'Mean Kmer Histogram Cluster {i}', fontsize=13)
        plt.savefig(f'mean-histogram_cluster_{i}.png')
        plt.close()

    # Plot clustering result in UMAP space
    fig, ax =plt.subplots()
    ax.scatter(clusterable_embedding[:,0], clusterable_embedding[:,1], c=self.clusterlabel, label=self.clusterlabel, cmap='Paired')
    ax.legend()
    ax.set_title('Normalized UMAP embedding of input reads (colored by HDBSCAN cluster)')
    plt.savefig('normalized_umap_clusterplot.png')
    plt.close()
    fig, ax = plt.subplots()
    ax = umap.plot.points(visual_embedding, labels=self.clusterlabel)
    ax.set_title(f'UMAP embedding colored by HDBSCAN cluster (ARI={adj_rs})')
    plt.savefig('umap_clusterplot.png')
    plt.close()

  def get_centroids(self, proc):
    """
    """
    seqCluster = { x : [] for x in set(self.clusterlabel)}
    for idx, cluster in self.allCluster:
      seqCluster[cluster].append(idx)

    p = Pool(self.proc)

    for cluster, sequences in seqCluster.items():
      if cluster == -1:
        continue 

      subProfiles = {seq : profile for seq, profile in Clusterer.d_profiles.items() if seq in sequences}

      for result in p.map(self.calc_pd, itertools.combinations(subProfiles.items(), 2)):
        seq1, seq2, dist = result
        Clusterer.matrix[seq1][seq2] = dist
        Clusterer.matrix[seq2][seq1] = dist

      tmpMinimum = math.inf
      centroidOfCluster = -1

      if len(sequences) == 1:
        centroidOfCluster = cluster[0]
        self.centroids.append(centroidOfCluster)
        continue

      for sequence in sequences:
        averagedDistance = 0

        for neighborSequence in sequences:
          if sequence == neighborSequence:
            continue
          averagedDistance += Clusterer.matrix[sequence][neighborSequence]

        averagedDistance /= len(sequences)-1

        if averagedDistance < tmpMinimum:
          tmpMinimum = averagedDistance
          centroidOfCluster = sequence

      self.centroids.append(centroidOfCluster)
    p.close()
    p.join()

  def output_centroids(self):
    """
    """
    reprSeqs = {centroid : self.d_sequences[centroid] for centroid in self.centroids}
    outputPath = f'{self.outdir}/{os.path.splitext(os.path.basename(self.sequenceFile))[0]}_hdbscan.fasta'

    with open(outputPath, 'w') as outStream:
      for centroidID, sequence in reprSeqs.items():
        outStream.write(f">{Clusterer.id2header[centroidID]}\n{sequence}\n")

    with open(f'{outputPath}.clstr', 'w') as outStream:
      for i in set(self.clusterlabel):
        outStream.write(f">Cluster {i}\n")
        seqInCluster = set([idx for idx,label in self.allCluster if label == i])
        for idx, seqIdx in enumerate(seqInCluster):
          outStream.write(f"{idx}\t{len(self.d_sequences[seqIdx])}nt, >{Clusterer.id2header[seqIdx]} ")
          if seqIdx in self.centroids:
            outStream.write('*\n')
          else:
            outStream.write("at +/13.37%\n")





########################################################
# CLUSTER PIPELINE
########################################################
def __abort_cluster(clusterObject, filename):
    logger.warn(f"Too few sequences for clustering in {os.path.basename(filename)}. No subcluster will be created.")
    del clusterObject

def perform_clustering():

  multiPool = Pool(processes=proc)
  virusClusterer = Clusterer(logger, inputSequences, lineageDict, primerDict, outdir,  k, proc, metric, neighbors, threshold, dimension, clusterSize, minSample)

  logger.info("Determining k-mer profiles for all sequences.")
  virusClusterer.determine_profile(multiPool)
  logger.info("Read lineage annotations for all sequences")
  virusClusterer.get_lineage_annotations()
  logger.info("Clustering with UMAP and HDBSCAN.")
  code = virusClusterer.apply_umap()
  #if code == 1:
  #  __abort_cluster(virusClusterer, inputSequences)
  #  return 0
  clusterInfo = virusClusterer.clusterlabel
  logger.info(f"Summarized {virusClusterer.dim} sequences into {clusterInfo.max()+1} clusters. Filtered {np.count_nonzero(clusterInfo == -1)} sequences due to uncertainty.")

  # logger.info("Extracting centroid sequences and writing results to file.\n")
  # logger.info("..get_centroids..\n")
  # virusClusterer.get_centroids(multiPool)
  # logger.info("..output_centroids..\n")
  # virusClusterer.output_centroids()

  logger.info(f"Extracting representative sequences for each cluster.")
  sequences = virusClusterer.d_sequences
  distanceMatrix = virusClusterer.matrix
  profiles = virusClusterer.d_profiles
  del virusClusterer

  return 0



#########################################################
# READ INPUT + LOGGER FUNCTIONALITIES
########################################################

inputSequences = None
outdir = None
k = None
proc = None
OUT = ''

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
    formatter = ColoredFormatter("%(log_color)sViralClust %(levelname)s -- %(asctime)s -- %(message)s", "%Y-%m-%d %H:%M:%S",
                                    log_colors={
                                            'DEBUG':    'bold_cyan',
                                            'INFO':     'bold_white',
                                            'WARNING':  'bold_yellow',
                                            'ERROR':    'bold_red',
                                            'CRITICAL': 'bold_red'}
                                )
    handle.setFormatter(formatter)
    logger.addHandler(handle)
    return logger

def create_outdir(outdir):
    try:
      os.makedirs(outdir)
      logger.info(f"Creating output directory: {outdir}")
    except FileExistsError:
      logger.warning(f"The output directory exists. Files will be overwritten.")

def parse_arguments(d_args):
  """
  Parse all given arguments and check for error (e.g. file names).

  Arguments:
  d_args -- dict with input parameters and their values

  Returns:
  parsed and checked parameter list
  """

  if d_args['--version']:
    print("viralClust version 0.1")
    exit(0)


  verbose = d_args['--verbose']
  if verbose:
    logger.setLevel(logging.INFO)

  inputSequences = d_args['<inputSequences>']
  print('inputSequences', inputSequences)
  if not os.path.isfile(inputSequences):
    logger.error("Couldn't find input sequences. Check your file")
    exit(1)

  lineageDict = d_args['<lineageDict>']
  print('lineageDict', lineageDict)
  if not os.path.isfile(lineageDict):
      logger.error("Couldn't find lineage dictionary. Check your file")
      exit(1)

  primerDict = d_args['<primerDict>']
  print('primerDict', primerDict)
  if not os.path.isfile(primerDict):
      logger.error("Couldn't find primer dictionary. Check your file")
      exit(1)

  try:
    k = int(d_args['--kmer'])
    print('k', k)
  except ValueError:
    logger.error("Invalid parameter for k-mer size. Please input a number.")
    exit(2)
    
  try:
    proc = int(d_args['--process'])
    print('proc', proc)
  except ValueError:
    logger.error("Invalid number for CPU cores. Please input a number.")
    exit(2)

  output = d_args['--output']
  print('output', output)
  if output == 'pwd':
    output = os.getcwd()
  now = str(datetime.now()).split('.')[0].replace(' ','_').replace(':','-')
  create_outdir(output)

  METRICES = [
              'euclidean', 'manhatten', 'chebyshev',
              'minkwoski', 'canberra', 'braycurtis',
              'mahalanobis', 'wminkowski',
              'seuclidean', 'cosine'
  ]
  metric = d_args['--metric']
  print('metric', metric)
  if metric not in METRICES:
    log.warn(f"Invalid metric chosen. Will default back to cosine distance.")
    metric = 'cosine'
  
  try:
    neighbors = int(d_args['--neighbors'])
    print('neighbors', neighbors)
    if neighbors <= 0:
      raise ValueError
      
  except ValueError:
      log.error("Invalid parameter for --neighbors. Please input a positive integer.")
      exit(2)

  try:
    threshold = float(d_args['--dThreshold'])
    print('threshold', threshold)
    if not 0.0 <= threshold < 1:
      raise ValueError
  except ValueError:
    log.error("Invalid parameter for --dThreshold. Please input a number between [0,1).")
    exit(2)

  try:
    dimension = int(d_args['--dimension'])
    print('dimension', dimension)
    if dimension < 1:
      raise ValueError
  except ValueError:
      log.error("Invalid parameter for --dimension. Please input a positive integer.")
      exit(2)
  
  try:
    clusterSize = int(d_args['--clusterSize'])
    print('clusterSize', clusterSize)
    if clusterSize < 1:
      raise ValueError
  except ValueError:
      log.error("Invalid parameter for --clusterSize. Please input a positive integer.")
      exit(2)

  if d_args['--minSample'] == "CLUSTERSIZE":
    minSample = clusterSize
    print('minSample', minSample)
  else:
    try:
      minSample = int(d_args['--minSample'])
      print('minSample', minSample)
      if minSample < 1:
        raise ValueError
    except ValueError:
        log.error("Invalid parameter for --minSample. Please input a positive integer.")
        exit(2)

  try:
    clusterThreshold = float(d_args['--clusterThreshold'])
    print('clusterThreshold', clusterThreshold)
    if clusterThreshold < 0:
      raise ValueError
  except ValueError:
      log.error("Invalid parameter for --clusterThreshold. Please input a positive integer.")
      exit(2)


  return (inputSequences, lineageDict,  primerDict, output,  k, proc, metric, neighbors, threshold, dimension, clusterSize, minSample, clusterThreshold)



###########################################################################
# MAIN
###########################################################################
if __name__ == "__main__":
  logger = create_logger()
  (inputSequences, lineageDict, primerDict, outdir, k, proc, metric, neighbors, threshold, dimension, clusterSize, minSample, clusterThreshold) = parse_arguments(docopt(__doc__))

  logger.info("Starting to cluster you data. Stay tuned.")
  perform_clustering()
  if os.path.islink(f"{os.path.dirname(outdir)}/latest"):
    os.remove(f"{os.path.dirname(outdir)}/latest")
  os.system(f"ln -s {outdir} {os.path.dirname(outdir)}/latest")

