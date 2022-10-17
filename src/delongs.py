#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare two ROC curves using DeLong's test (DeLong, E.R. et al. 1988)

@Returns z-statistic and p-value from DeLong's test

@author: marionfischer
"""

import numpy, pandas
import scipy.stats as st
from argparse import ArgumentParser
from sklearn import metrics

"""
# Create the command line parser
parser = ArgumentParser(
    This script performs the DeLong's test (DeLong, E.M. et al. 1988)
    to compare two ROC curves and find the z-score and p-score
    indicating if the two curves are similar (p>0.05) or not (p<=0.05).
)
parser.add_argument(
    "-a",
    "--method_A_predictions",
    dest="predsA",
    help=INPUT file containing array of size n of n predicted values
    using method A
)
parser.add_argument(
    "-b",
    "--method_B_predictions",
    dest="predsB",
    help=INPUT file containing array of size n of n predicted values
    using method B
)
parser.add_argument(
    "-c",
    "--classifier_labels",
    dest="actual",
    help=INPUT file containing array of size n of n expected or
    classifier labels (e.g. 0 or 1)
)

parser.add_argument(
    "-s",
    "--summary",
    dest="summary",
    help=Output file to append summary statistics
)

# Parse the command line
args = parser.parse_args()
"""
# Calculate AUC (= Mann-Whitney statistic applied to samples X and Y)
def auc(X, Y):
    auc = 0
    if len(X)==0 or len(Y)==0:
        auc = 0.5
    else:
        auc = 1/(len(X)*len(Y)) * sum([kernel(x,y) for x in X for y in Y])
    return auc

# Implements step function to compare two instances of opposing classes
def kernel(X, Y):
    return 0.5 if Y==X else int(Y < X)

# Estimate variance (V) and covariance (C)
def structural_components(X, Y):
    V10 = [1/len(Y) * sum([kernel(x,y) for y in Y]) for x in X]
    V01 = [1/len(X) * sum([kernel(x,y) for x in X]) for y in Y]
    return V10, V01

# Generates an estimated variance-covariance matrix (S)
def get_S_entry(V_A, V_B, auc_A, auc_B):
    return 1/(len(V_A)-1) * sum([(a-auc_A)*(b-auc_B) for a,b in zip(V_A,V_B)])

# Calculates z-score
def z_score(var_A, var_B, covar_AB, auc_A, auc_B):
    return (auc_A - auc_B)/((var_A + var_B - 2*covar_AB)**(0.5))

# Separate prediction vectors by classifier label
def group_by_label(preds, actual):
    X = [p for (p,a) in zip(preds,actual) if a]
    Y = [p for (p,a) in zip(preds,actual) if not a]
    return X, Y
# Compute entries of covariance matrix S (covar_AB=covar_BA)
def computeCoVarMatrix(V_A10, V_A01, V_B10, V_B01, auc_A, auc_B):
    if len(V_A10)==0 or len(V_A01)==0 or len(V_B10)==0 or len(V_B01)==0:
        var_A = numpy.nan
        var_B = numpy.nan
        covar_AB = numpy.nan
    else:
        var_A = (get_S_entry(V_A10, V_A01, auc_A, auc_A) * 1/len(V_A10)
             +get_S_entry(V_A01, V_A01, auc_A, auc_A) * 1/len(V_A01))
        var_B = (get_S_entry(V_B10, V_B01, auc_B, auc_B) * 1/len(V_B10)
             +get_S_entry(V_B01, V_B01, auc_B, auc_B) * 1/len(V_B01))
        covar_AB = (get_S_entry(V_A10, V_B10, auc_A, auc_B) * 1/len(V_A10)
                +get_S_entry(V_A01, V_B01, auc_A, auc_B) * 1/len(V_A01))
    return var_A, var_B, covar_AB

def main():

    #__________Read input arrays__________
    #preds_A = numpy.genfromtxt(args.predsA)
    #preds_B = numpy.genfromtxt(args.predsB)

    #__________Compute variance and AUC by label__________
    print("Computing variance and auc")
    X_A, Y_A = group_by_label(preds_A, actual)
    X_B, Y_B = group_by_label(preds_B, actual)

    V_A10, V_A01 = structural_components(X_A, Y_A)
    V_B10, V_B01 = structural_components(X_B, Y_B)

    auc_A = auc(X_A, Y_A)
    auc_B = auc(X_B, Y_B)

    #__________Compute entries of covariance matrix S (covar_AB=covar_BA)__________
    print("Computing covariance matrix")
    var_A, var_B, covar_AB = computeCoVarMatrix(V_A10, V_A01, V_B10, V_B01, auc_A, auc_B) 

    #__________Two tailed test__________
    if numpy.isnan(var_A) is False or numpy.isnan(var_B) is False or numpy.isnan(covar_AB) is False:
        z = z_score(var_A, var_B, covar_AB, auc_A, auc_B)
        p = st.norm.sf(abs(z))*2
    else:
        z = numpy.nan
        p = 1

    #__________Print output__________
    #with open(args.summary, 'a') as s:
    #    s.write(str(z)+" "+str(p)+"\n")
    print("Done")
if __name__ == "__main__":
    main()
