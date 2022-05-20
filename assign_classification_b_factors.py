#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: marionfischer
"""

from argparse import ArgumentParser
import numpy
import pandas
import re

# Create the command line parser
parser = ArgumentParser(
    """
    This script returns multiple files containing the input need for 
    insert_b_factor.py for each pdb where the threshold identified using
    Youden's J statistic provides the cutoff for assigning if a prediction
    is negative (below threshold) or positive (above threshold) and if the 
    prediction is true or false
    """
)
parser.add_argument(
    "-r",
    "--ROC_curve",
    dest="roc",
    help="""Input file containing ROC curve data to calculate
    Youden's J statistic with first 3 columns:
    <fpr> <tpr> <threshold> <virus>
    """
)

parser.add_argument(
    "-p",
    "--pdb_list",
    dest="pdblist",
    help="""Input file containing list of pdb_bfactor output file names and #residues per pdb in columns:
    <PDB_b_factor.txt> <#residues of PDB>
    Note: list must be in same order as PDB output for NeighborSum#.txt
    """
)
parser.add_argument(
    "-j",
    "--j_statistic",
    dest="summary",
    required=False,
    help="""Output file to write (append) J statistic"""
)
parser.add_argument(
    "-v",
    "--prediction_values",
    dest="values",
    help="""Input file (NeighborSum#.txt) containing prediction/predictor values"""
)

# Parse the command line
args = parser.parse_args()

def calculate_Youdens_J( rocs ):
    """
    @Returns Youden's J statistic from given ROC curves and separates by
    test set (left out viral protein ensemble)
    """
    
    #----------Read in ROC curve data----------
    ROC = pandas.read_csv(rocs, delimiter=' ', header=None)
    fpr = ROC.iloc[:,0].to_numpy()       # false positive rate (1 - specificity)
    tpr = ROC.iloc[:,1].to_numpy()       # true positive rate (sensitivity)
    thresholds = ROC.iloc[:,2].to_numpy() # threshold values
    viruses = ROC.iloc[:,3].to_list()
    
    #----------Calculate J statistic-----------
    # J = max(sensitivity + specificity - 1)
    J = 0.0 
    J_threshold = 0.0
    for v in range( len(tpr) ):
        j_ind = tpr[v] + (1-fpr[v]) - 1
        if j_ind > J:
            J = j_ind
            J_threshold = thresholds[v]
    print ("J:",J, " threshold:",J_threshold)
    return J, J_threshold

def assign_classification_label( J_threshold, prediction, expected ):
    """
    J_threshold = Youden's J statistic given ROC data
    prediction_data = NeighborSum#.txt (file containing calculated prediction and predictor values
    column = column to use for assigning labels within NeighborSum#.txt (prediction values used to calculate ROC)

    @returns array of size n (with n=len(prediction_data) of labels 
    based on if prediction is a true positive(0), 
    true negative(3), false positive(1), or false negative(2)
    """

    labels = numpy.zeros(len(prediction))
    for i in range( len(prediction)):
        # True positives
        if prediction[i] >= J_threshold and expected[i] == 1.0:
            labels[i] = int(0)
        # False positives
        if prediction[i] >= J_threshold and expected[i] == 0.0:
            labels[i] = int(1)
        # False negatives
        if prediction[i] < J_threshold and expected[i] == 1.0:
            labels[i] = int(2)
        # True negatives
        if prediction[i] < J_threshold and expected[i] == 0.0:
            labels[i] = int(3)
            
    return labels

def  create_b_factor_resdata( pdb ):
    """
    @ returns array of size (n,3) containing necessary input for inserting the classification labels
    as B factors using the 'insert_b_factor.py script with the following format:
    <chain> <res#> <label>
    with n being the number of residues in a given pdb model
    """

    r = []
    c = []
    # Parse PDB
    pdb_parsed = []
    with open(pdb, 'r') as p:
        pdb = p.readlines()
        for line in pdb:
            # Account for numerous instances where white space no longer exists between columns
            # Residue numbers >=1000
            num_index = 22
            num_line = line[:num_index] + ' ' + line[num_index:]
            # Y coordinate (add 1 index from initial index 22)
            y_index = 39
            y_line = num_line[:y_index] + ' ' + num_line[y_index:]
            # Z coordinate (add 2 index from initial index 46)
            z_index = 48
            z_line = y_line[:z_index] + ' ' + y_line[z_index:]
            pdb_parsed.append(re.split('\s+', z_line))
    for res in enumerate(pdb_parsed):
        if res[1][0]=="ATOM" and res[1][2]=="CA":
            c.append(res[1][4])
            r.append(res[1][5])
            
    return c, r

def main():

    #---------------Calculate J statistic-----------------------------
    # Lists of J and threshold values for each virus
    J, J_threshold = calculate_Youdens_J( args.roc )

    #---------------Create B factor files for each pdb----------------
    pdb_data = pandas.read_csv(args.pdblist, header=None, delimiter = ' ')
    pdb_names = pdb_data.iloc[:, 0].to_list()
    pdb_output = pdb_data.iloc[:, 1].to_list()
    pdb_lengths = pdb_data.iloc[:,3].to_numpy()
    prediction_values = numpy.genfromtxt(args.values, usecols=(1))
    expected_values = numpy.genfromtxt(args.values, usecols=(0))

    # Iterate through file containing all predictor and expected values
    predictor_index = 0
    for p in range( len(pdb_names) ):
        print("Assigning labels for", pdb_names[p])
        res_predictors = numpy.zeros(pdb_lengths[p])
        res_expected = numpy.zeros(pdb_lengths[p])
        for r in range( pdb_lengths[p]):
            res_predictors[r] = prediction_values[predictor_index]
            res_expected[r] = expected_values[predictor_index]
            predictor_index += 1
            pdb_labels = assign_classification_label( J_threshold, res_predictors, res_expected)
        chain, res_num = create_b_factor_resdata( pdb_names[p] )
        with open(pdb_output[p], 'w') as b:
                 for l in range( pdb_lengths[p]):
                     b.write(str(chain[l])+" "+str(res_num[l])+" "+str(pdb_labels[l])+"\n")

    print("Done.")

# Invoke main
main()
