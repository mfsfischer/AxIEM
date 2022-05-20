#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Linear regression model for epitope prediction using predictor values only

@Returns output file containing precision and recall curve
@Returns output file containig ROC curve
@Returns output file containing summary statistics

@author: marionsauer
"""

import math, numpy, pandas
from argparse import ArgumentParser
from sklearn.linear_model import LinearRegression
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from Bio.PDB.PDBParser import PDBParser
from collections import Counter

# Create the command line parser
parser = ArgumentParser(
    """This script returns a per-residue score based on the sum of all
    predictor scores of neighboring residues and itself as well as the 
    precision/recall values."""
)
parser.add_argument(
    "-a",
    "--training-dataset",
    dest="training",
    help="""INPUT space-separated file containing training dataset.
    Requires second to last column to contain PDB file (path) name and 
    last column to contain output target variable as a one hot encoded
    variable"""
)
parser.add_argument(
    "-b",
    "--testing-dataset",
    dest="testing",
    help="""INPUT space-separated file containing training dataset.
    Requires second to last column to contain PDB file (path) name and 
    last column to contain output target variable as a one hot encoded
    variable"""
)
parser.add_argument(
    "-e",
    "--epitope_radius",
    dest="radius",
    help="""Radius of maximum distance away from reference residue that
    could form an epitope"""
)
parser.add_argument(
    "-l",
    "--epitope_lower",
    dest="lower",
    help="""Radius of minimum distance away from reference residue that
    could form an epitope with a weight of 1"""
)
parser.add_argument(
    "-p",
    "--precision_recall",
    dest="prec_recall",
    help="""Output file to append precision recall curve"""
)
parser.add_argument(
    "-r",
    "--roc_curves",
    dest="roc",
    help="""Output file to append precision recall curve"""
)
parser.add_argument(
    "-z",
    "--summary",
    dest="summary",
    help="""Output file to append summary statistics"""
)
parser.add_argument(
    "-y",
    "--predicted_values",
    dest="results",
    help="""Output file to append predicted values of epitopes"""
)

# Parse the command line
args = parser.parse_args()
# Parse PDB file
pdbparser = PDBParser()

def readDataset( data ):
    """

    Parameters
    ----------
    input_data : file containing per-residue predictors (pred#x), PDB file 
    name,and per-residue classifier, with the following format:
        <pred#1> <pred#2> <...pred#n> </path/to/name>.pdb> <classifier[0/1]>

    Returns
    -------
    A list of predictors of size (res,p), where res is the number of residues
    in the dataset and p is the number predictors assigned to each res, a list
    of pdb names of size res, and a list containing classifier labels of size (res)

    """
    X_in = pandas.DataFrame(data.iloc[:, :-2]).to_numpy()
    pdbs = data.iloc[:,-2].tolist()
    y_in = data.iloc[:,-1].to_numpy()
    
    # -----Extract residue predictors by PDB------
    pdb_counts = dict(Counter(pdbs).items())    # Get number, lengths, names of PDBs in dataset
    num_pdbs = len(pdb_counts)                  # Number of PDBs
    pdb_lengths = list(pdb_counts.values())     # Lengths of PDBs
    pdb_names = list(pdb_counts.keys())         # Names of PDBs
    
    pdb_index = 0
    res_index = 0
    X = numpy.zeros((X_in.shape[0], X_in.shape[1]+1))
    y = numpy.zeros(y_in.shape[0])
    for p in range(num_pdbs):
        pdb_data = []
        pdb = pdb_names[p]
        classifier = []
        for res in range(pdb_lengths[p]):
            pdb_data.append(X_in[pdb_index])
            classifier.append(y_in[pdb_index])
            pdb_index += 1
        print("Setting (non-)epitopes for", pdb, "with", len(pdb_data), "residues")
        structure = pdbparser.get_structure("apdb", pdb)
        pdb_neighbor_predictors = getResidueScores(structure, numpy.asarray(pdb_data))
        for res in range(pdb_neighbor_predictors.shape[0]):
            X[res_index] = pdb_neighbor_predictors[res]
            y[res_index] = classifier[res]
            res_index += 1
    
    return X, y

def getResidueScores(structure, predictors):
    """
    Find score of each residue as linear sum of predictors. Recalculate  
    residue score as the weighted sum of linear scores where weights correspond to 
    relative distance. Returns predictor set with neighbor-weighted score
    as an additional predictor.
    """
    # Find linear sum of all predictors
    scores = numpy.zeros( predictors.shape[0] )
    for res in range( 0, predictors.shape[0] ):
        for p in range( 0, predictors.shape[1] ):
            scores[res] += predictors[res][p]
    
    # Reweight score to account for distance to each residue
    # and return new predictors as 
    # {p1, p2, ..., pn, neighbor_score}
    X = numpy.zeros( (predictors.shape[0], predictors.shape[1]+1) )
    for res in range( 0, len(scores) ):
        # Get weighted neighborhood weights for reference residue
        wts = getNeighborWeights( structure, res )
        neighbor_score = 0
        for w in range(0, len(wts) ):
            neighbor_score += wts[w] * scores[w]
        for p in range( predictors.shape[1]+1 ):
            if p < predictors.shape[1]:
                X[res][p] = predictors[res][p]
            if p == predictors.shape[1]:
                X[res][p] = neighbor_score
    return X

def getNeighborWeights(structure, residue_num ):
    """
    Defines weight a cosine function of distance to reference residue
    Any residues within 4 Angstroms of reference residue's Cbeta atom is considered
    direct neighbors and is given a weight of 1. Any greater than 4 and at
    less than r Angstroms from ref Cbeta atom is weighted sigmoidally.
    Any greater distance that is not considered a neighbor and given a weight
    of 0. Returns weights of all other residue distances in relation to the
    input residue ID number
    """
    r = float(args.radius)
    l = int(args.lower)
    residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
    wts = numpy.empty( len(residues) )
    for x in range(0, len(residues)):
        cb_distance = 0.0
        # Check for glycine for reference residue
        # If glycine, use virtual C-beta atom coordinate
        if residues[residue_num].get_resname()=='GLY':
            cb1 = numpy.array( getVirtualCbeta( residues[residue_num] ) )
        else:
            cb1 = numpy.array( residues[residue_num]["CB"].get_coord() )
        # Check for glycine for comparison residue 
        if residues[x].get_resname()=='GLY':
            cb2 = numpy.array( getVirtualCbeta( residues[x] ) )
        else:
            cb2 = numpy.array( residues[x]["CB"].get_coord() )
        cb_distance = numpy.linalg.norm(cb2 - cb1)
        if cb_distance <= l:
            wts[x] = 1.0
        if cb_distance > l and cb_distance <= r:
            distance_bound_ratio  =  ( ( cb_distance - l ) / (r - l ) ) * math.pi
            sigmoidal_proximity = 0.5 * ( math.cos( distance_bound_ratio ) + 1 )
            wts[x] = sigmoidal_proximity
        else:
            wts[x] = 0.0
    return wts

def getVirtualCbeta( glycine ):
    """
    Gets a glycine's N, C, and CA coordinates as vectors to find a rotation
    matrix that rotates N -120 degrees along the CA-C vector, with the origin
    being the CA coordinates
    """
    n = glycine["N"].get_coord()
    c = glycine["C"].get_coord()
    ca = glycine["CA"].get_coord()
    # Place origin at C-alpha atom
    n = n - ca
    c = c - ca
    # Find rotation matrix that rotates n -120 degrees along the 
    # ca-c vector
    theta = -120.0/180.0 * math.pi
    rot = rotaxis( theta, c )
    # Apply rotation to ca-n vector
    cb_origin = numpy.dot( n, rot )
    return cb_origin + ca

def rotaxis( theta , vec ):
    """Calculate left multiplying rotation matrix
    Taken from Bio.PDB.vectors.rotaxis2m since rotaxis2m vector was being
    assigned as a numpy.ndarray vectors object not as a rotaxis2m parameter

    Parameters
    ----------
    theta : type float. The rotation angle.
    vec : type array. The rotation axis.

    Returns
    -------
    rot : The rotation matrix, a 3x3 array.

    """
    vec = normalize(vec)
    c = numpy.cos(theta)
    s = numpy.sin(theta)
    t = 1 - c
    x, y, z = numpy.array( vec )
    rot = numpy.zeros((3,3))
    # 1st row
    rot[0, 0] = t * x * x + c
    rot[0, 1] = t * x * y - s * z
    rot[0, 2] = t * x * z + s * y
    # 2nd row
    rot[1, 0] = t * x * y + s * z
    rot[1, 1] = t * y * y + c
    rot[1, 2] = t * y * z - s * x
    # 3rd row
    rot[2, 0] = t * x * z - s * y
    rot[2, 1] = t * y * z + s * x
    rot[2, 2] = t * z * z + c
    return rot

def normalize( vec ):
    norm = numpy.sqrt(sum(vec * vec))
    vec_normalized = vec / norm
    return vec_normalized

def main():
    # -----------------Input Data-------------------
    training_data = pandas.read_csv(args.training, delimiter=' ', header=None)
    testing_data = pandas.read_csv(args.testing, delimiter=' ', header=None)
    
    # -----------------Extract Data-------------------
    print("Getting data with training set", args.training,"...")
    train_X, train_y = readDataset( training_data )
    test_X, test_y = readDataset( testing_data )
    
    # -----------------Train Model-------------------
    model = LinearRegression(normalize=True).fit(train_X, train_y)
    r_sq = model.score(train_X, train_y)
    print('coefficient of determination:', r_sq)
    print('intercept:', model.intercept_)
    print('slope:', model.coef_)
    
    # -----------------Predict with Model and Evaluate-------------------
    prediction = model.predict(test_X)
    precision, recall, pr_threshold = precision_recall_curve(test_y, prediction)
    fpr, tpr, roc_threshold = roc_curve(test_y, prediction)
    auc_score = roc_auc_score(test_y, prediction)
    print("AUC score:", auc_score)
    
    # -----------------Write Output-------------------
    with open( args.prec_recall, 'a') as pr, open( args.roc, 'a') as r, open( args.results, 'a') as result, open(args.summary, 'a') as z:
        for i in range( len(pr_threshold) ):
            pr.write( str(precision[i])+" "+str(recall[i])+" "+str(pr_threshold[i])+" "+args.testing+" NeighborSum "+args.radius+"\n")
        for j in range( len(roc_threshold) ):
            r.write( str(fpr[j])+" "+str(tpr[j])+" "+str(roc_threshold[j])+" "+args.testing+" NeighborSum "+args.radius+"\n")
        for k in range( len(test_y) ):
            result.write(str(test_y[k])+" "+str(prediction[k])+" "+str(test_X[k][0])+" "+str(test_X[k][1])+" "+str(test_X[k][2])+" "+str(test_X[k][3])+" "+args.testing+" NeighborSum "+args.radius+"\n")
        z.write(str(r_sq)+" "+str(model.intercept_)+" "+str(model.coef_)+" "+str(auc_score)+" "+args.training+" "+args.testing+" NeighborSum "+args.radius+"\n")

# Invoke main()
main()
