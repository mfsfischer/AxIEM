#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 10:16:12 2020

@author: marionsauer
"""
from argparse import ArgumentParser
import numpy, pandas
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

# Create the command line parser
parser = ArgumentParser(
    """Returns the precision-recall and ROC curves, as well as the AUC value
    given a 1D array of predicted values and a 1D array, both of size n,
    provided as two input files."""
)
parser.add_argument(
    "-d",
    "--data",
    dest="data",
    help="""NeighborSum#.txt from output running LinearRegression_NeighborSum.py 
    """
)
parser.add_argument(
    "-e",
    "--radius",
    dest="radius",
    help="""Radius used to compute predictors
    """
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

# Parse the command line
args = parser.parse_args()

def main():
    # -----------------Input Data-------------------
    #dataset = pandas.read_csv(args.data, header=None, delimiter=' ')
    #expected = dataset.iloc[:,0].to_numpy()
    #predicted = dataset.iloc[:,1].to_numpy()
    expected = numpy.genfromtxt(args.data, delimiter=' ', usecols=(0))
    predicted = numpy.genfromtxt(args.data, delimiter=' ', usecols=(1))
    
    
    # -----------------Evaluate Data-------------------
    precision, recall, pr_threshold = precision_recall_curve(expected, predicted)
    fpr, tpr, roc_threshold = roc_curve(expected, predicted)
    auc_score = roc_auc_score(expected, predicted)
    print(args.radius, auc_score)
    
    # -----------------Write Output-------------------
    with open( args.prec_recall, 'a') as pr, open( args.roc, 'a') as r:
        for i in range( len(pr_threshold) ):
            pr.write( str(precision[i])+" "+str(recall[i])+" "+str(pr_threshold[i])+"\n")
        for j in range( len(roc_threshold) ):
            r.write( str(fpr[j])+" "+str(tpr[j])+" "+str(roc_threshold[j])+" "+ args.data+"\n")

# Invoke main()
main()
