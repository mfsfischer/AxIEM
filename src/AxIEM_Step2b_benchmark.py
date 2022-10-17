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
from sklearn.preprocessing import MinMaxScaler

# Create command line parser
parser = ArgumentParser(
    """Step 2 in AxIEM benchmark:
    Normalizes data using the MinMaxScaler to bounds [0,1]
    """
)
parser.add_argument(
    "--data",
    dest="data",
    required=True,
    help="""INPUT from Step 1's output features text file - see documentation for formatting"""
)
parser.add_argument(
    "--scaled",
    dest="scaled",
    help="""OUTPUT file to write all scaled feature values to be used in Step 3- see documentation for formatting"""
)

# Parse the command line
args = parser.parse_args()

def scaleData(data):
    """
    Parameters
    ----------
    data        : text file of updated features from Step 2 - see documentation

    Returns
    -------
    scaled_data : numpy array size (n, 6) of scaled features
    viruses     : list size (n) of virus name labels
    pdbs        : list size (n) of PDB name labels
    resno       : list size (n) of residue number ids
    classifiers : list size (n) of classifier labels
    """
    data = pandas.read_csv(data, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    pdbs = data.iloc[:,1].to_list()
    classifiers = data.iloc[:,2].to_numpy()
    resno = data.iloc[:,3].to_numpy()
    X = pandas.DataFrame(data.iloc[:,4:-1]).to_numpy()
    scaler = MinMaxScaler()
    scaled_data = scaler.fit_transform(X)
    return scaled_data, viruses, pdbs, classifiers, resno

def main():

    print("\nStep 2b: Scale features to [0,1]")
    scaled_features, viruses, pdbs, classifiers, resno = scaleData(args.data)

    # Write output to use as scaled features for testing models in Step 3
    with open(args.scaled, 'w') as update:
        update.write("Virus PDB Classifier ResNo REU CP NV reu8 cp8 nv8 reu16 cp16 nv16 reu24 cp24 nv24 reu32 cp32 nv32 reu40 cp40 nv40 reu48 cp48 nv48 reu56 cp56 nv56 reu64 cp64 nv64 reuR cpR nvR\n")
        for i in range(scaled_features.shape[0]):
            update.write('%r %r %r %r ' % (viruses[i], pdbs[i], classifiers[i], resno[i]))
            for j in range(scaled_features.shape[1]):
                if j==scaled_features.shape[1]-1:
                    update.write(str(scaled_features[i][j])+"\n")
                else:
                    update.write(str(scaled_features[i][j])+" ")

    print("Done")

if __name__ == "__main__":
    main()
