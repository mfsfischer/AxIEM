#!/usr/bin/env python3
"""
Gets AUC value given expected classifier label and prediction score
@author marionfischer
"""
# Basic calculations and statistical analysis
from argparse import ArgumentParser
import numpy
import pandas
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

# Create command line parser
parser = ArgumentParser(
    """Gets AUC value given expected classifier label column and prediction score column"""
    )
parser.add_argument(
    "--data",
    dest="data",
    required=True,
    help="""INPUT text file; must be space delimited"""
    )
parser.add_argument(
    "--classifier",
    dest="classifier",
    required=True,
    help="""Column # with index=0 of classifier labels""")
parser.add_argument(
    "--predictions",
    dest="prediction",
    required=True,
    help="""Column # with index=0 of prediction scores""")

# Parse the command line
args = parser.parse_args()

def main():
    data = pandas.read_csv(args.data, delimiter=' ')
    classifier = pandas.DataFrame(data.iloc[:,int(args.classifier)])
    prediction = pandas.DataFrame(data.iloc[:,int(args.prediction)])
    fpr, tpr, threshold = roc_curve(classifier, prediction, pos_label=1)
    auc = roc_auc_score(classifier, prediction)
    print("AUC value:", auc)

if __name__=='__main__':
    main()
