#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:38:53 2020

@author: marionfischer
"""

from argparse import ArgumentParser
import numpy
import math

# Create the command line parser
parser = ArgumentParser(
    """
    This script returns the Youden's J statistic of a given ROC
    curve to identify the threshold that optimizes sensitivity
    and specificity
    """
)
parser.add_argument(
    "-d",
    "--data",
    dest="data",
    help="""File containing expected and predicted values as columns <expected> <predicted>"""
)
parser.add_argument(
    "-r",
    "--ROC_curve",
    dest="roc",
    help="""Input file containing ROC curve data with first 3 columns:
    <fpr> <tpr> <threshold>
    """
)
parser.add_argument(
    "-s",
    "--statistic_output",
    dest="summary",
    help="""Output file to write (append) J and MCC statistics"""
)

# Parse the command line
args = parser.parse_args()

def main():

    #----------Read in ROC curve data----------
    ROC = numpy.genfromtxt(args.roc, delimiter=' ', usecols=(0,1,2))
    fpr = numpy.array(ROC[:,0])       # false positive rate (1 - specificity)
    tpr = numpy.array(ROC[:,1])       # true positive rate (sensitivity)
    threshold = numpy.array(ROC[:,2]) # threshold values
    
    #----------Calculate J statistic-----------
    # J = max(sensitivity + specificity - 1)
    J = 0.0
    J_threshold = 0.0
    
    for i in range(0,tpr.shape[0]):
        j_ind = tpr[i] + (1-fpr[i]) - 1
        if j_ind > J:
            J = j_ind
            J_threshold = threshold[i]

    print("J: ", J)

    #----------Calculate MCC statistic---------
    expected = numpy.genfromtxt(args.data, delimiter=' ', usecols=(0))
    predicted = numpy.genfromtxt(args.data, delimiter=' ', usecols=(1))
    TP = 0
    FP = 0
    TN = 0
    FN = 0
    MCC = 0
    for i in range(len(expected)):
        if predicted[i] >= J_threshold and expected[i] == 1.0:
            TP += 1
        if predicted[i] >= J_threshold and expected[i] == 0.0:
            FP += 1
        if predicted[i] < J_threshold and expected[i] == 1.0:
            FN += 1
        if predicted[i] < J_threshold and expected[i] == 0.0:
            TN += 1
    numerator = (TP*TN) - (FP*FN)
    if (TP+FP)==0 or (TP+FN)==0 or (TN+FP)==0 or (TN+FN)==0:
        MCC = 0
        if (TP+FP)==0:
            print("TP+FP==0")
        elif (TP+FN)==0:
            print("TP+FN==0")
        elif (TN+FP)==0:
            print("TN+FP==0")
        elif (TN+FN)==0:
            print("TN+FN==0")
    else:
        denominator = math.sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
        MCC = numerator / denominator
    print("MCC: ", MCC)
    #----------Write output--------------------
    with open(args.summary, 'a') as s:
        s.write(args.roc + " " + str(J) + " " + str(J_threshold) + " " + str(MCC) +"\n")

# Invoke main()
main()
