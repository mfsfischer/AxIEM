#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Linear regression model for epitope prediction using predictor values only

@Returns output file containing precision and recall curve
@Returns output file containig ROC curve
@Returns output file containing summary statistics

@author: marionsauer
"""

import numpy, pandas
from argparse import ArgumentParser
from sklearn.linear_model import LinearRegression
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score


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
parser.add_argument(
    "-e",
    "--epitope_radius",
    dest="radius",
    help="""Upper Boundary used for epitope radius size"""
)

# Parse the command line
args = parser.parse_args()

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
    X = data.iloc[:, 5].to_numpy()
    X = X.reshape(-1,1)
    y = data.iloc[:, 0].to_numpy()
    return X, y

def main():
    # -----------------Input Data-------------------
    training_data = pandas.read_csv(args.training, delimiter=' ', header=None)
    testing_data = pandas.read_csv(args.testing, delimiter=' ', header=None)
    
    # -----------------Extract Data-------------------
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
            pr.write( str(precision[i])+" "+str(recall[i])+" "+str(pr_threshold[i])+" "+str(args.testing)+" NS "+args.radius+"\n")
        for j in range( len(roc_threshold) ):
            r.write( str(fpr[j])+" "+str(tpr[j])+" "+str(roc_threshold[j])+" "+str(args.testing)+" NS "+args.radius+"\n")
        for k in range( len(test_y) ):
            result.write(str(test_y[k])+" "+str(prediction[k])+" "+str(args.testing)+" NS "+args.radius+"\n")
        z.write(str(r_sq)+" "+str(model.intercept_)+" "+str(model.coef_)+" "+str(auc_score)+" "+str(args.testing)+" NS "+args.radius+"\n")

# Invoke main()
main()
