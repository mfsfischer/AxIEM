#!/usr/bin/env python3
"""
This script identifies residues that may be "mislabeled" as an epitope for a given conformation as well as identifies residues predicted to be epitopes between aligned proteins. 
"""

from argparse import ArgumentParser
import numpy
import pandas
import sys

parser = ArgumentParser()
parser.add_argument(
    "--protein1",
    dest="protein1",
    required=True,
    help="""Input space-delimited file containing matrix of size n by c+1 matrix with n=#aligned residues and c=#protein conformations. File should formatted as such: <Residue #> <label 1> <label2> ... <label_c> ;Note, if alignment contains gaps, gaps should be indicated with 'NaN' """
)
parser.add_argument(
    "--protein2",
    dest="protein2",
    required=False,
    help="""Optional: Input space-delimited file containing matrix of size n by c matrix with n=#aligned residues and c=#protein conformations, where n=#residues and n(protein1)=n(protein2). Format should be the same as for --protein1 input"""
)
parser.add_argument(
    "--output",
    dest="output",
    required=True,
    default='epitope_findings.txt',
    help="""Output file name. Default is 'epitope_findings.txt'"""
)

args = parser.parse_args()

def findPossibleMislabels( labels ):
    """
    @params LABELS : array of residue classification labels of size (n,c) with n =  number of residues and c = number of conformations per protein
    @return incorrect fraction false positives (nFP)
    @return incorrect fraction false negatives (nFN)
    @return array of size n of classifications, with 0=not an epitope in all cases or gap, 1=epitope in at least 1 case
    """

    FP_count = 0
    FN_count = 0
    nFP = 0        # False positives when also false negatives
    nFN1 = 0       # False negatives when also true positives 
    nFN2 = 0       # False negatives when also false positives
    classifications = numpy.zeros(labels.shape[0])
    residues = 0

    for n in enumerate(labels):
        # Iterate to count FP and FN
        # If for a residue there is either a TP&FN, FP&FN then classification as an epitope was mis-assigned to the wrong conformation
        # Not counting TN&FP, as TN cannot be defined
        FP = 0
        FN = 0
        TP = 0
        for c in enumerate(n[1]):
            if c[1]==0.0 or c[1]==1.0 or c[1]==2.0:
                if c[1]==0.0:
                    TP +=1
                elif c[1]==1.0:
                    FP+=1
                    FP_count+=1
                elif c[1]==2.0:
                    FN+=1
                    FN_count+=1
            else:
                # Assumes that an indicated gap within a protein of aligned conformation sequences contains same number of gaps
                classifications[n[0]]=int(0)
        # Count number of times a FN was predicted when in another conformation the same residue was predicted as a TP,
        # i.e. that conformation should also have been assigned as an expected epitope
        if TP>0:
            classifications[n[0]]=int(1)
            if FN>0 and FP==0:
                nFN1+=FN
        # Count number of times expected outcomes were switched (FP and FN > 0); Count as epitope residue
        elif FN > 0 and FP > 0:
            nFP+=FP
            nFN2+=FN
            classifications[n[0]]=int(1)
    if FP_count != 0:
        nFP = nFP/FP_count
    else:
        nFP = 0
    if FN_count != 0:
        nFN1 = nFN1/FN_count
        nFN2 = nFN2/FN_count
    else:
        nFN1 = 0
        nFN2 = 0
    print("Number of False Positives:",FP_count)
    print("Fraction inverted:",nFP)
    print("Number of False Negatives",FN_count)
    print("Fraction incorrect:",nFN1)
    print("Fraction inverted:",nFN2,"\n")
    return nFP, nFN1, nFN2, classifications

def findAlignedEpitopes( classifications1, classifications2, res1, res2 ):
    """
    @params classifications1 : array of size n containing binary classification labels from findPossibleMislabels for protein1
    @params classifications2 : array of size n containing binary classification labels ^ for aligned protein2
    @params res1 : array of size (n,2) containing residue#, chainID for protein1
    @params res2 : array of size (n,2) containing residue#, chainID for protein2
    @returns epitopes : list of tuples containing [ (protein1 res#, chainID), (protein2 res#, chainID)] of positions where both are classified as epitopes (1)
    """
    if len(classifications1) != len(classifications2):
        print("Proteins to find aligned epitopes are not of same alignment length.")
        print("Protein 1 has)",len(classifications1),"positions while protein 2 has",len(classifications2),"positions")
        print("Aborting...")
        sys.exit(-1)
    else:
        # Get residue number for aligned epitope residues
        epitopes = []
        for position in range( len(classifications1)):
            if classifications1[position]==1 and classifications2[position]==1:
                epitopes.append((res1[position], res2[position]))
            else:
                continue
    return epitopes

def main():

    # Parse input data
    protein1 = pandas.read_csv(args.protein1, delimiter=' ', header=None)
    protein2 = pandas.read_csv(args.protein2, delimiter=' ', header=None)
    labels1 = protein1.iloc[:, 2:].to_numpy()
    labels2 = protein2.iloc[:, 2:].to_numpy()
    res1 = protein1.iloc[:, 0:2].values.tolist()
    res2 = protein2.iloc[:, 0:2].values.tolist()

    # Calculate percentage possible incorrect FP,FN and update classification labels
    nfp_1, nfn1_1, nfn2_1, classifications1 = findPossibleMislabels(labels1)
    nfp_2, nfn1_2, nfn2_2, classifications2 = findPossibleMislabels(labels2)

    # Find aligned epitope residue ids in protein1 and protein2
    epitopes = findAlignedEpitopes(classifications1, classifications2, res1, res2)

    # Write output
    with open(args.output, 'w') as o:
        o.write("Protein FN_TP FP_switchedFN FN_switchedFP\n")
        o.write("1 "+str(nfn1_1)+" "+str(nfp_1)+" "+str(nfn2_1)+"\n")
        o.write("2 "+str(nfn1_2)+" "+str(nfp_2)+" "+str(nfn2_2)+"\n")
        o.write("\n\nAligned epitope residues\n")
        for e in epitopes:
            o.write(str(e)+"\n")
            
main()
