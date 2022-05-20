#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Linear regression model for epitope prediction using predictor values only
plus randomly chosen value for Neighbor Sum

@Returns output file containing precision and recall curve
@Returns output file containig ROC curve
@Returns output file containing summary statistics

@author: marionsauer
"""

import math, numpy, random
from argparse import ArgumentParser
import pomegranate
from pomegranate import *
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

# Create the command line parser
parser = ArgumentParser(
    """This script returns a per-residue randomized score based on the sum of all
    predictor scores of neighboring residues and itself as well as the 
    precision/recall values."""
)
parser.add_argument(
    "-a",
    "--training-dataset",
    dest="training",
    help="""INPUT space-separated file containing training dataset.
    (NeighborSumValues#.txt of training set)"""
)
parser.add_argument(
    "-b",
    "--testing-dataset",
    dest="testing",
    help="""INPUT space-separated file containing training dataset.
    (NeighborSumValues#.txt of test set)"""
)
parser.add_argument(
    "-e",
    "--epitope_radius",
    dest="radius",
    help="""Radius of maximum distance away from reference residue that
    could form an epitope"""
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

def randomizeDataset( data ):
    """

    Parameters
    ----------
    input_data : file containing per-residue predictors (pred#x), PDB file 
    name with the following format:
        <classifier> <predictionvalue> <pred#1> <pred#2> <...pred#n> </path/to/name>.pdb> <> <radius>

    Returns
    -------
    A list of predictors of size (res,p), where res is the number of residues
    in the dataset and p is the number predictors assigned to each res, a list
    of pdb names of size res, and a list containing classifier labels of size (res)

    """
    X_in = numpy.genfromtxt(data, delimiter=' ', usecols=(2,3,4,5))
    y = numpy.genfromtxt(data, delimiter=' ', usecols=(0))
    
    # Find mean and standard deviation of each predictor's values
    mu = numpy.mean(X_in, axis=1)
    sigma = numpy.std(X_in, axis=1)

    # Assign random gaussian value for Neighbor Sum based on input mu, sigma
    # where mu and sigma are from mean and std dev of predictor distributions
    # from algorithm using same boundary conditions and feature sets
    X = numpy.zeros((X_in.shape[0], X_in.shape[1]))
    for residue in range( X_in.shape[0]):
        for predictor in range( X_in.shape[1]):
            X[residue][predictor] = random.gauss(mu[predictor], sigma[predictor])
    return X, y

def main():
    # -----------------Randomize Data-------------------
    train_X, train_y = randomizeDataset( args.training )
    test_X, test_y = randomizeDataset( args.testing )

    # -----------------Train Model-------------------
    # With initial + NS
    model = pomegranate.BayesClassifier.from_samples(pomegranate.distributions.MultivariateGaussianDistribution, train_X, train_y)
    print("Bayes classifier mean prediction: ", (model.predict(test_X) == test_y).mean())
    
    # -----------------Predict with Model and Evaluate-------------------
    model_predict = model.predict(test_X)
    model_probabilities = model.predict_proba(test_X)
    prediction = model_probabilities[:,1]
    precision, recall, pr_threshold = precision_recall_curve(test_y, prediction)
    fpr, tpr, roc_threshold = roc_curve(test_y, prediction)
    auc_score = roc_auc_score(test_y, prediction)
    print("AUC score:", auc_score)
    
    # -----------------Write Output-------------------
    with open( args.prec_recall, 'a') as pr, open(args.summary, 'a') as z:
        for i in range( len(pr_threshold) ):
            pr.write( str(precision[i])+" "+str(recall[i])+" "+str(pr_threshold[i])+" "+args.training+" random "+args.radius+"\n")
    with open( args.roc, 'a') as r:       
        for i in range( len(roc_threshold) ):
            r.write( str(fpr[i])+" "+str(tpr[i])+" "+str(roc_threshold[i])+" "+args.training+" random "+args.radius+"\n")
    with open( args.results, 'a') as result:
        for k in range( len(test_y) ):
            result.write(str(test_y[k])+" "+str(prediction[k])+" "+str(test_X[k][0])+" "+str(test_X[k][1])+" "+str(test_X[k][2])+" "+str(test_X[k][3])+" "+args.training+" random "+args.radius+"\n")
    with open( args.summary, 'a') as z:    
        z.write(args.training+" Random "+args.radius+" "+str(auc_score)+"\n")
# Invoke main()
main()
