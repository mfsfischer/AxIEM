#!/usr/bin/env python3
"""
Performs leave-one-out analyses of epitope mapping models using a range of epitope radii
@author marionfischer
"""
# Parsing files
from argparse import ArgumentParser

# Basic calculations and statistical analysis
from collections import Counter, defaultdict
import math
import numpy
import pandas
import random

# Create command line parser
parser = ArgumentParser(
    """Step 2 in AxIEM benchmark:
    Computes 'neighborhood' feature scores for multiple upper boundary conditions and
    outputs computed feature values into file (to avoid having to compute values for
    each leave-one-out test)
    """
)
parser.add_argument(
    "--data",
    dest="data",
    required=True,
    help="""INPUT from Step 1's output features text file - see documentation for formatting"""
)
parser.add_argument(
    "--features",
    dest="features",
    help="""OUTPUT file to write all feature values to be used in Step 3- see documentation for formatting"""
)
parser.add_argument(
    "--randomized_features",
    dest="random",
    help="""OUTPUT file to write all randomized feature values to be used in Step 3"""
)
parser.add_argument(
    "--plotting",
    dest="plotting",
    help="""OUTPUT file to write all feature value formatted for plotting distributions"""
)

# Parse the command line
args = parser.parse_args()

def getNeighborSumWts_byVirus(dataset):
    """
    Defines weight a cosine function of distance to reference residue
    Any residues within 4 Angstroms of reference residue's Cbeta atom is considered
    direct neighbors and is given a weight of 1. Any greater than 4 and at
    less than u Angstroms from reference Cbeta atom is weighted sigmoidally.
    Any greater distance that is not considered a neighbor and given a weight
    of 0.

    Parameters
    ----------
    dataset : txt file of annotated data (e.g. AxIEM.data)

    Returns
    -------
    viruses_wts : dict of dict containing numpy arrays size (9,n,n) of each virus's PDB's Neighbor Sum by upper boundary weights
    viruses : list len(n) of virus names - used for writing output
    pdbs : list len(n) of PDB IDs - used for writing output
    classifiers : list len(n) of classifier labels - used for writing output
    resno : list len(n) of residue IDs - used for writing output
    """
    data = pandas.read_csv(dataset, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    pdbs = data.iloc[:,1].to_list()
    classifiers = data.iloc[:,2].to_list()
    resno = data.iloc[:,3].to_list()
    viruses_wts = defaultdict(lambda : defaultdict([]))
    for i in range(len(viruses)):
        if not viruses[i] in viruses_wts:
            viruses_wts[viruses[i]] = {}
            if not pdbs[i] in viruses_wts[viruses[i]]:
                cb_distance_map = numpy.genfromtxt("pdb_contact_maps/{}maps".format(pdbs[i]))
                pdblength = len(cb_distance_map)
                viruses_wts[viruses[i]][pdbs[i]] = getNSwts(cb_distance_map, pdblength)
        else:
            if not pdbs[i] in viruses_wts[viruses[i]]:
                cb_distance_map = numpy.genfromtxt("pdb_contact_maps/{}maps".format(pdbs[i]))
                viruses_wts[viruses[i]][pdbs[i]] = getNSwts(cb_distance_map, pdblength)
                pdblength = len(cb_distance_map)
    return viruses_wts, viruses, pdbs, classifiers, resno

def getNSwts(cb_distances, pdblength):
    """
    Parameters
    ----------
    cb_distances : numpy array size (n,n) of single PDB's contact map distances
    pdblength : int, number of residues in protein to determine 'R' boundary

    Returns
    -------
    wts : numpy array size (9,n,n) of single PDB's wts given Neighbor Sum's upper boundary limit
    """
    upper = ['8','16','24','32','40','48','56','64','R']
    l = 4.0
    wts = numpy.empty((9,cb_distances.shape[0], cb_distances.shape[1]))
    for x in range(cb_distances.shape[0]):
        for y in range(cb_distances.shape[1]):
            # Like Cb distances, wts are symmetric
            # Iterate u values in each case to minimize iterating cb_distances
            if y == x:
                for u in range(len(upper)):
                    wts[u][x][y] = 1.0
            if y > x:
                if cb_distances[x][y] <= l:
                    for u in range(len(upper)):
                        wts[u][x][y] = 1.0
                        wts[u][y][x] = 1.0
                else:
                    i = 0
                    for u in upper:
                        u = getUpperValue(u, pdblength)
                        if cb_distances[x][y] > l and cb_distances[x][y] <= u:
                            distance_bound_ratio  =  ( ( cb_distances[x][y] - l ) / (u - l ) ) * math.pi
                            sigmoidal_proximity = 0.5 * ( math.cos( distance_bound_ratio ) + 1 )
                            wts[i][x][y] = sigmoidal_proximity
                            wts[i][y][x] = sigmoidal_proximity
                        else:
                            wts[i][x][y] = 0.0
                            wts[i][y][x] = 0.0
                        i += 1
    return wts

def getUpperValue(u, pdblength):
    """
    Parameters
    ----------
    u : str, where 'R' is the pdb length dependent upper boundary radius
    pdblength : int
    
    Returns
    -------
    u : float
    """
    if u != 'R':
        u = float(u)
    else:
        u = 15.1 + 0.0147 * pdblength
    return u

def getPDBfeatures(dataset, virus, pdb):
    """
    Parameters
    ----------
    dataset : args.data

    Returns
    -------
    pdb_features : numpy array size (n,3) of a single PDB's per-residue features
    """
    data = pandas.read_csv(dataset, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    pdbs = data.iloc[:,1].to_list()
    features = data.iloc[:,4:7].to_numpy()
    pdb_features = []
    for i in range(len(viruses)):
        if viruses[i]==virus and pdbs[i]==pdb:
            #pdb_features = numpy.concatenate((pdb_features, features[i]))
            pdb_features.append(features[i])
    return numpy.asarray(pdb_features)

def insertNeighborSums(wts, features):
    """
    Calculate Neighborhood Sum scores for each residue's per-residue features.

    Parameters
    ----------
    wts : numpy array size (u,n,n) of single PDB's NS_u weights
    features : numpy array size (n,f) with n=number of residues within a protein and f={REU,CP,NV}

    Returns
    -------
    updated_features : numpy array size (n, f+f*u) of original features + Neighbor Sums of each feature for each tested upper boundary radius
    """
    updated_features = numpy.zeros( (features.shape[0], features.shape[1]+features.shape[1]*wts.shape[0]) )
    for res in range(features.shape[0]):
        # For each upper boundary radius, calculate Neighbor Sum value for each per-residue feature
        for u in range(wts.shape[0]):
            neighbor_score_reu = 0
            neighbor_score_cp  = 0
            neighbor_score_nv  = 0
            for wt in range(wts.shape[1]):
                neighbor_score_reu += wts[u][res][wt] * features[res][0]
                neighbor_score_cp  += wts[u][res][wt] * features[res][1]
                neighbor_score_nv  += wts[u][res][wt] * features[res][2]
            updated_features[res][features.shape[1]+u*features.shape[1]+0] = neighbor_score_reu
            updated_features[res][features.shape[1]+u*features.shape[1]+1] = neighbor_score_cp
            updated_features[res][features.shape[1]+u*features.shape[1]+2] = neighbor_score_nv
        for p in range(features.shape[1]):
            updated_features[res][p] = features[res][p]
    return numpy.asarray(updated_features)

def randomize( features ):
    """
    Parameters
    ----------
    features : numpy array size (n,f) of input feature value array

    Returns
    -------
    randomized_features = numpy array size (n,f) of features that have been randomized using the gaussian distribution's mu and sigma of each feature
    """
    randomized_features = numpy.zeros((features.shape[0], features.shape[1]))
    mu = numpy.mean(features, axis=0)
    sigma = numpy.std(features, axis=0)
    for residue in range(features.shape[0]):
        for feature in range(features.shape[1]):
            randomized_features[residue][feature] = random.gauss(mu[feature], sigma[feature])
    return randomized_features

def main():

    print("\nStep 2: Add Neighbor Sum values to data set \n 1)Calculating Neighbor Sum weights")

    # Need each PDB's NSu weights and feature sets grouped by viral protein for NS calculations
    # where NSu is the upper boundary radius (u) dependent Neighbor Sum (NS) that is used to
    # calculate each residue's cumulative neighbor-weighted REU, CPvar, and NV scores
    NSu_PDBwts_byVirus, viruses, pdbs, classifiers, resno = getNeighborSumWts_byVirus(args.data) 

    # For each viral protein's PDB, calculate all NSu-weighted features to insert NSu scores to input dataset
    print("2) Adding Neighbor Sum weighted features")
    updated_features = []
    for virus, pdb in NSu_PDBwts_byVirus.items():
        for p, ns_wts in pdb.items():
            pdb_features = getPDBfeatures(args.data, virus, p)
            pdb_updated_features = insertNeighborSums(ns_wts, pdb_features)
            for i in range(len(pdb_updated_features)):
                updated_features.append(pdb_updated_features[i])
    updated_features = numpy.asarray(updated_features)
    randomized_features = randomize(updated_features)

    print("3) Writing output")
    # Write output to use as features for testing models in Step 3
    with open(args.features, 'w') as update:
        update.write("Virus PDB Classifier ResNo REU CP NV reu8 cp8 nv8 reu16 cp16 nv16 reu24 cp24 nv24 reu32 cp32 nv32 reu40 cp40 nv40 reu48 cp48 nv48 reu56 cp56 nv56 reu64 cp64 nv64 reuR cpR nvR\n")
        for i in range(updated_features.shape[0]):
            update.write('%r %r %r %r ' % (viruses[i], pdbs[i], classifiers[i], resno[i]))
            for j in range(updated_features.shape[1]):
                if j==updated_features.shape[1]-1:
                    update.write(str(updated_features[i][j])+"\n")
                else:
                    update.write(str(updated_features[i][j])+" ")
    # Write output to use as randomized features for negative control in testing models in Step 3
    with open(args.random, 'w') as rand:
        rand.write("Virus PDB Classifier ResNo REU CP NV reu8 cp8 nv8 reu16 cp16 nv16 reu24 cp24 nv24 reu32 cp32 nv32 reu40 cp40 nv40 reu48 cp48 nv48 reu56 cp56 nv56 reu64 cp64 nv64 reuR cpR nvR\n")
        for i in range(randomized_features.shape[0]):
            rand.write('%r %r %r %r ' % (viruses[i], pdbs[i], classifiers[i], resno[i]))
            for j in range(randomized_features.shape[1]):
                if j==randomized_features.shape[1]-1:
                    rand.write(str(randomized_features[i][j])+"\n")
                else:
                    rand.write(str(randomized_features[i][j])+" ")
    # Write output for plotting... I know it's clunky     
    with open(args.plotting, 'w') as plot:
        plot.write("Virus PDB Classifier REU CP NV Uweight")
        u = 0
        for i in range(updated_features.shape[0]):
            line = [viruses[i], pdbs[i], classifiers[i], updated_features[i][u+0], updated_features[i][u+1], updated_features[i][u+2], "Per-residue"]
            plot.write(" ".join(map(str, line)))
            plot.write("\n")
        u = 3
        for i in range(updated_features.shape[0]):
            line = [viruses[i], pdbs[i], classifiers[i], updated_features[i][u+0], updated_features[i][u+1], updated_features[i][u+2], "8"]
            plot.write(" ".join(map(str, line)))
            plot.write("\n")
        u = 6
        for i in range(updated_features.shape[0]):
            line = [viruses[i], pdbs[i], classifiers[i], updated_features[i][u+0], updated_features[i][u+1], updated_features[i][u+2], "16"]
            plot.write(" ".join(map(str, line)))
            plot.write("\n")
        u = 9
        for i in range(updated_features.shape[0]):
            line = [viruses[i], pdbs[i], classifiers[i], updated_features[i][u+0], updated_features[i][u+1], updated_features[i][u+2], "24"]
            plot.write(" ".join(map(str, line)))
            plot.write("\n")
        u = 12
        for i in range(updated_features.shape[0]):
            line = [viruses[i], pdbs[i], classifiers[i], updated_features[i][u+0], updated_features[i][u+1], updated_features[i][u+2], "32"]
            plot.write(" ".join(map(str, line)))
            plot.write("\n")
        u = 15
        for i in range(updated_features.shape[0]):
            line = [viruses[i], pdbs[i], classifiers[i], updated_features[i][u+0], updated_features[i][u+1], updated_features[i][u+2], "40"]
            plot.write(" ".join(map(str, line)))
            plot.write("\n")
        u = 18
        for i in range(updated_features.shape[0]):
            line = [viruses[i], pdbs[i], classifiers[i], updated_features[i][u+0], updated_features[i][u+1], updated_features[i][u+2], "48"]
            plot.write(" ".join(map(str, line)))
            plot.write("\n")
        u = 21
        for i in range(updated_features.shape[0]):
            line = [viruses[i], pdbs[i], classifiers[i], updated_features[i][u+0], updated_features[i][u+1], updated_features[i][u+2], "56"]
            plot.write(" ".join(map(str, line)))
            plot.write("\n")
        u = 24
        for i in range(updated_features.shape[0]):
            line = [viruses[i], pdbs[i], classifiers[i], updated_features[i][u+0], updated_features[i][u+1], updated_features[i][u+2], "64"]
            plot.write(" ".join(map(str, line)))
            plot.write("\n")
        u = 27
        for i in range(updated_features.shape[0]):
            line = [viruses[i], pdbs[i], classifiers[i], updated_features[i][u+0], updated_features[i][u+1], updated_features[i][u+2], "R"]
            plot.write(" ".join(map(str, line)))
            plot.write("\n")
    print("Done")

if __name__ == "__main__":
    main()
