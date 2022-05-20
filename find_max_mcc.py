#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: marionfischer
"""

from argparse import ArgumentParser
import numpy
import math

# Create the command line parser
parser = ArgumentParser(
    """
    This script identifies the maximum Matthew's correlation coefficient
    and the associated threshold given an ROC curve and prediction scores
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
    help="""Output file to write (append) associated threshold and max MCC statistics"""
)

# Parse the command line
args = parser.parse_args()

def getMCC(threshold):
    """
    @params threshold : score threshold to define predicted positives and negatives
    @return MCC : Matthew's correlation coefficient calculated given the threshold
    """
    expected = numpy.genfromtxt(args.data, delimiter=' ', usecols=(0))
    predicted = numpy.genfromtxt(args.data, delimiter=' ', usecols=(1))
    TP = 0
    FP = 0
    TN = 0
    FN = 0
    MCC = 0
    for i in range(len(expected)):
        if predicted[i] >= threshold and expected[i] == 1.0:
            TP += 1
        if predicted[i] >= threshold and expected[i] == 0.0:
            FP += 1
        if predicted[i] < threshold and expected[i] == 1.0:
            FN += 1
        if predicted[i] < threshold and expected[i] == 0.0:
            TN += 1
    numerator = (TP*TN) - (FP*FN)
    if (TP+FP)==0 or (TP+FN)==0 or (TN+FP)==0 or (TN+FN)==0:
        MCC = 0
    else:
        denominator = math.sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
        MCC = numerator / denominator
    return MCC
    

def main():

    #----------Read in ROC curve data----------------------------
    ROC = numpy.genfromtxt(args.roc, delimiter=' ', usecols=(0,1,2))
    fpr = numpy.array(ROC[:,0])       # false positive rate (1 - specificity)
    tpr = numpy.array(ROC[:,1])       # true positive rate (sensitivity)
    threshold = numpy.array(ROC[:,2]) # threshold values

    #----------Get max MCC statistic and associated threshold----
    threshold_mcc = 0
    MCC = 0
    for t in enumerate(threshold):
        mcc = getMCC(t[1])
        if mcc > MCC:
            MCC = mcc
            print("updated mcc", mcc)
            threshold_mcc = t[1]

    #----------Write output--------------------------------------
    print("MCC: ", MCC)
    print("threshold: ", threshold_mcc)
    with open(args.summary, 'a') as s:
        s.write(args.roc + str(threshold_mcc) + " " + str(MCC) +"\n")

# Invoke main()
main()
