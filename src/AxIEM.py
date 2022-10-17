#!/usr/bin/env python3
"""
Performs leave-one-out analyses of epitope mapping models using a range of epitope radii
@author marionfischer
"""
# Parsing files
from argparse import ArgumentParser
from Bio.PDB.PDBParser import PDBParser
from Bio import BiopythonWarning

# Basic calculations and statistical analysis
from collections import Counter, defaultdict
import math
import numpy
import pandas
import pomegranate
from pomegranate import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
import random

# Compare performance ROC curves
import scipy.stats as st
import delongs

# Mute Biopython PDB parsing warnings when reading Rosetta output files
import warnings
warnings.simplefilter('ignore', BiopythonWarning)

# Create command line parser
parser = ArgumentParser(
    """
    Runs AxIEM features with option to use linear regression, Bayes classifier, or both
    """
    )
parser.add_argument(
    "--data",
    dest="data",
    required=True,
    help="""INPUT text file containing pre-computed REU, CPrmsd, and NV per-residue features
    that shares same format as AxIEM.data - see documentation on how to generate file"""
)
parser.add_argument(
    "--features",
    dest="feature_sets",
    help="""OUTPUT file to write residue feature set values including Neighbor Sum scaled feature sets"""
)
parser.add_argument(
    "--virus_label",
    dest="label",
    help="""virus label as string"""
)
parser.add_argument(
    "--model",
    dest="model",
    help="""Model to use for prediction; options include 'LR' for linear regression, or 'BC' for Bayes classifier"""
)

# Parse the command line
args = parser.parse_args()
# Parse PDB file
pdbparser = PDBParser()

def getPDBdistances_byVirus(dataset, test_protein):
    """
    Parameters
    ----------
    dataset : txt file of annotated data (e.g. AxIEM.data)

    Returns
    -------
    viruses_cb_distances : dict of dict containing numpy arrays size (n,n) of each virus's PDB's contact map distances
    R : float, upper boundary radius to use to calculate Neighbor Sum
    """
    data = pandas.read_csv(dataset, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    pdbs = data.iloc[:,1].to_list()
    R = 0.0
    viruses_cb_distances = defaultdict(lambda : defaultdict([]))
    for i in range(len(viruses)):
        if not viruses[i] in viruses_cb_distances:
            viruses_cb_distances[viruses[i]] = {}
            if not pdbs[i] in viruses_cb_distances[viruses[i]]:
                pdb = "pdb_structures/{}".format(pdbs[i])
                structure = pdbparser.get_structure("apdb", pdb)
                residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
                print("Getting {} structure {} with {} residues".format(viruses[i], pdbs[i], len(residues)))
                if viruses[i]==test_protein:
                    R = 15.1 + 0.0147 * len(residues)
                viruses_cb_distances[viruses[i]][pdbs[i]] = getCbDistances(residues)
        else:
            if not pdbs[i] in viruses_cb_distances[viruses[i]]:
                pdb = "pdb_structures/{}".format(pdbs[i])
                structure = pdbparser.get_structure("apdb", pdb)
                residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
                print("Getting {} structure {} with {} residues".format(viruses[i], pdbs[i], len(residues)))
                viruses_cb_distances[viruses[i]][pdbs[i]] = getCbDistances(residues)
    return viruses_cb_distances, R

def getCbDistances(residues):
    """
    Parameters
    ----------
    residues : list of len(n) of Biopython PDB structure's residue entities

    Returns
    -------
    cb_distances : numpy array size (n,n) of all C-beta atoms' pairwise distances
    """
    cb_distances = numpy.zeros((len(residues), len(residues)))
    for res in range(len(residues)):
        cb1 = []
        if residues[res].get_resname()=='GLY':
            cb1 = numpy.array( getVirtualCbeta(residues[res]) )
        else:
            cb1 = numpy.array( residues[res]["CB"].get_coord() )
        for d in range(len(residues)):
            # Cb contact map is symmetric - only compute half of matrix
            cb2 = []
            if d > res:
                if residues[d].get_resname()=='GLY':
                    cb2 = numpy.array( getVirtualCbeta(residues[d]) )
                else:
                    cb2 = numpy.array( residues[d]["CB"].get_coord() )
                cb_distances[res][d] = numpy.linalg.norm(cb2 - cb1)
                cb_distances[d][res] = numpy.linalg.norm(cb2 - cb1)
            elif d == res:
                cb_distances[res][d] = 0.0       
    return cb_distances

def getVirtualCbeta( glycine ):
    """
    Gets a glycine's N, C, and CA coordinates as vectors to find a rotation
    matrix that rotates N -120 degrees along the CA-C vector, with the origin
    being the CA coordinates
    """
    n = glycine["N"].get_coord()
    c = glycine["C"].get_coord()
    ca = glycine["CA"].get_coord()
    # Place origin at C-alpha atom
    n = n - ca
    c = c - ca
    # Find rotation matrix that rotates n -120 degrees along the 
    # ca-c vector
    theta = -120.0/180.0 * math.pi
    rot = rotaxis( theta, c )
    # Apply rotation to ca-n vector
    cb_origin = numpy.dot( n, rot )
    return cb_origin + ca

def rotaxis( theta , vec ):
    """Calculate left multiplying rotation matrix
    Taken from Bio.PDB.vectors.rotaxis2m since rotaxis2m vector was being
    assigned as a numpy.ndarray vectors object not as a rotaxis2m parameter

    Parameters
    ----------
    theta : type float. The rotation angle.
    vec : type array. The rotation axis.

    Returns
    -------
    rot : The rotation matrix, a 3x3 array.

    """
    vec = normalize(vec)
    c = numpy.cos(theta)
    s = numpy.sin(theta)
    t = 1 - c
    x, y, z = numpy.array( vec )
    rot = numpy.zeros((3,3))
    # 1st row
    rot[0, 0] = t * x * x + c
    rot[0, 1] = t * x * y - s * z
    rot[0, 2] = t * x * z + s * y
    # 2nd row
    rot[1, 0] = t * x * y + s * z
    rot[1, 1] = t * y * y + c
    rot[1, 2] = t * y * z - s * x
    # 3rd row
    rot[2, 0] = t * x * z - s * y
    rot[2, 1] = t * y * z + s * x
    rot[2, 2] = t * z * z + c
    return rot

def normalize( vec ):
    norm = numpy.sqrt(sum(vec * vec))
    vec_normalized = vec / norm
    return vec_normalized

def getNeighborWeights_byVirus(wts, cb_distances, dataset, upper, test_protein):
    """
    Defines weight a cosine function of distance to reference residue
    Any residues within 4 Angstroms of reference residue's Cbeta atom is considered
    direct neighbors and is given a weight of 1. Any greater than 4 and at
    less than u Angstroms from reference Cbeta atom is weighted sigmoidally.
    Any greater distance that is not considered a neighbor and given a weight
    of 0.

    Parameters
    ----------
    wts : initialized empty list
    cb_distances : dict of dict containing numpy arrays size (n,n) of each virus's PDB's contact map distances
    dataset : txt file of annotated data (e.g. AxIEM.data)
    upper : int, value to set upper boundary for neighbor sum calculation

    Returns
    -------
    viruses_wts : dict of dict containing numpy arrays size (n,n) of each virus's PDB's Neighbor Sum weights 
    """
    data = pandas.read_csv(dataset, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    pdbs = data.iloc[:,1].to_list()
    viruses_wts = defaultdict(lambda : defaultdict([]))
    for i in range(len(viruses)):
        if not viruses[i] in viruses_wts:
            viruses_wts[viruses[i]] = {}
            if not pdbs[i] in viruses_wts[viruses[i]]:
                viruses_wts[viruses[i]][pdbs[i]] = getNSwts(cb_distances[viruses[i]][pdbs[i]], upper)
        else:
            if not pdbs[i] in viruses_wts[viruses[i]]:
                viruses_wts[viruses[i]][pdbs[i]] = getNSwts(cb_distances[viruses[i]][pdbs[i]], upper)
    return viruses_wts

def getNSwts(cb_distances, upper):
    """
    Parameters
    ----------
    cb_distances : numpy array size (n,n) of single PDB's contact map distances
    upper : int, value to set upper boundary for neighbor sum calculation

    Returns
    -------
    wts : numpy array size (n,n) of single PDB's wts given Neighbor Sum's upper boundary limit
    """
    u = upper
    l = 4.0
    wts = numpy.empty((cb_distances.shape[0], cb_distances.shape[1]))
    for x in range(cb_distances.shape[0]):
        for y in range(cb_distances.shape[1]):
            # Like Cb distances, wts are symmetric
            if y == x:
                wts[x][y] = 1.0
            if y > x:
                if cb_distances[x][y] <= l:
                    wts[x][y] = 1.0
                    wts[y][x] = 1.0
                if cb_distances[x][y] > l and cb_distances[x][y] <= u:
                    distance_bound_ratio  =  ( ( cb_distances[x][y] - l ) / (u - l ) ) * math.pi
                    sigmoidal_proximity = 0.5 * ( math.cos( distance_bound_ratio ) + 1 )
                    wts[x][y] = sigmoidal_proximity
                    wts[y][x] = sigmoidal_proximity
                else:
                    wts[x][y] = 0.0
                    wts[y][x] = 0.0
    return wts

def constructDatasets(data, virus_cb_distances, virus_wts, test_protein):
    """
    Parameters
    ----------
    data               : text file of annotated complete dataset with each line [VIRUS PDB CLASSIFIER_LABEL RES_ID REU CPrmsd NV] of len(n)
    virus_cb_distances : dict of dict containing numpy arrays size (n,n) of each virus's PDB's contact map distances
    virus_wts          : dict of dict containing numpy arrays size (n,n) of each virus's PDB's NS_u weights
    test_protein       : string of protein name to be left out of training set and to be used for testing

    Returns
    -------
    X         : numpy array list size (len(n-m), 4) of scaled training feature set + NS_u calculated with scaled data
    y         : numpy array size len(m) of classifier labels
    test_X    : numpy array list size (len(m), 4) of scaled test feature set + NS_u calculated with scaled data
    test_y    : numpy array size len(m) of classifier labels
    """
    data = pandas.read_csv(data, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    pdbs = data.iloc[:,1].to_list()
    allX = pandas.DataFrame(data.iloc[:,4:7]).to_numpy()
    scaler = MinMaxScaler()
    allXscaled = scaler.fit_transform(allX)
    allY = data.iloc[:,2].to_numpy()

    X = []
    y = []
    test_X = []
    test_y = []

    # Add Neighbor Sum (NS_u) feature using scaled per-residue features
    allX = assignNeighborSum( allXscaled, pdbs, viruses, virus_cb_distances, virus_wts)
    
    # Separate into training and test sets
    for v in range(len(viruses)):
        if viruses[v] == test_protein:
            test_X.append(allX[v])
            test_y.append(allY[v])
        else:
            X.append(allX[v])
            y.append(allY[v])
    return numpy.asarray(X), numpy.asarray(y), numpy.asarray(test_X), numpy.asarray(test_y)

def assignNeighborSum( X_in, pdbs, viruses, CBdicts, WTSdicts ):
    """
    Parameters
    ----------
    X_in : numpy array size (n, scaled) where scaled==MinMaxScaler([REU, CPrmsd, NV], axis=0)
    pdbs : list size (n) of pdb names
    viruses : list size (n) of virus names
    CBdicts : dict of dicts to get numpy array size (n,n) of cb_distances for each PDB
    upper : float value to use as upper boundary for NS_u feature calculated in getNeighborWeights

    Returns
    -------
    X : numpy array size (n, p) where p==[scaled-REU, scaled-CPrmsd, scaled-NV, NS_u with scaled features], axis=0)

    """
    # -----Extract residue predictors by PDB------
    pdb_counts = dict(Counter(pdbs).items())         # Get number, lengths, names of PDBs in dataset - all ordered
    num_pdbs = len(pdb_counts)                       # Number of PDBs
    pdb_lengths = list(pdb_counts.values())          # Lengths of PDBs
    pdb_names = list(pdb_counts.keys())              # Names of PDBs
    
    pdb_index = 0
    res_index = 0
    X = numpy.zeros((X_in.shape[0], X_in.shape[1]+1))
    # For each PDB in the dataset
    for p in range(num_pdbs):
        pdb_data = []
        pdb = "pdb_structures/{}".format(pdb_names[p])
        virus = viruses[res_index]
        pdb_cb_dists = CBdicts[virus][pdb_names[p]]
        pdb_wts = WTSdicts[virus][pdb_names[p]]
        # Get all residues within single PDB 
        for res in range(pdb_lengths[p]):
            pdb_data.append(X_in[pdb_index])
            pdb_index += 1
        # Add Neighbor Sum feature given weights calculated with a given upper boundary
        pdb_NS = addNeighborSum(pdb_cb_dists, pdb_wts, numpy.asarray(pdb_data))
        # Reconstruct entire data set with same indexing as X_in
        for res in range(pdb_NS.shape[0]):
            X[res_index] = pdb_NS[res]
            res_index += 1
    return X

def addNeighborSum(cb_distances, wts, features):
    """
    Add Neighborhood Sum score of each residue as weighted sum of raw/scaled feature sets.

    Parameters
    ----------
    cb_distances : numpy array size (n, n) of single PDB's contact map distances
    wts : numpy array size (n, n) of single PDB's NS_u weights
    features : numpy array size (n, f) with n=number of residues within a protein and f=number of features
    upper : float value to use as upper boundary for NS_u feature calculated in getNeighborWeights

    Returns
    -------
    X : numpy array size (n, f+1) 
    """
    # Find inner sum of all scaled per-residue features
    inner_sum = numpy.zeros(features.shape[0])
    for res in range(features.shape[0]):
        for p in range(features.shape[1]):
            inner_sum[res] += features[res][p]
    # Find outer sum of all surrounding weighted residues' features
    # accounting for all distances to each residue
    # and return new features as 
    # {f1, f2, ..., fn, neighbor_score} or
    X = numpy.zeros( (features.shape[0], features.shape[1]+1) )
    for res in range(len(inner_sum)):
        # Get weighted neighborhood weights for reference residue
        neighbor_score = 0
        for w in range(len(wts)):
            neighbor_score += wts[res][w] * inner_sum[res]
        for p in range(features.shape[1]+1):
            if p < features.shape[1]:
                X[res][p] = features[res][p]
            if p == features.shape[1]:
                X[res][p] = neighbor_score
    return X

# The following three functions include lines for model performance output (##)
# They have been silenced for brevity
def performLinearRegression( train_X, train_y, test_X, test_y):
    # -----------------Train Model-------------------
    # Line below - scikit-learn linear regression v 1.0 will soon be deprecated
    # model = LinearRegression(normalize=True).fit(train_X, train_y)
    # model is scaled but due to bimodal distributions of features, model shouldn't be normalized
    #model = make_pipeline(StandardScaler(with_mean=False), LinearRegression())
    model = LinearRegression().fit(train_X, train_y)
    ## r_sq = model.score(train_X, train_y)
    ## print('coefficient of determination:', r_sq)
    ## print('intercept:', model.intercept_)
    ## print('linear regression coefficients:', model.named_steps['linearregression'].coef_)
    ## coefficients = model.named_steps['linearregression'].coef_
    
    # -----------------Predict with Model and Evaluate-------------------
    prediction = model.predict(test_X)
    fpr, tpr, roc_threshold = roc_curve(test_y, prediction)
    auc = roc_auc_score(test_y, prediction)
    #print("Linear regression AUC:", auc)
    return auc, fpr, tpr, roc_threshold, prediction

def performBayesClassification(train_X, train_y, test_X, test_y):
    model = pomegranate.BayesClassifier.from_samples(pomegranate.distributions.MultivariateGaussianDistribution, train_X, train_y)
    ## print("Bayes classifier mean prediction: ", (Bayes_model.predict(test_X) ==  test_y ).mean())
    # Predict with Bayes model and evaluate
    predict = model.predict(test_X)
    probabilities = model.predict_proba(test_X)
    prediction = probabilities[:,1]
    fpr, tpr, roc_threshold = roc_curve(test_y, prediction)
    auc = roc_auc_score(test_y, prediction)
    #print("Bayes Classification AUC:", auc)
    return auc, fpr, tpr, roc_threshold, prediction

def performRandomForestClassification(train_X, train_y, test_X, test_y):
    model = RandomForestClassifier()
    model.fit(train_X, train_y)
    predict = model.predict(test_X)
    probabilities = model.predict_proba(test_X)
    prediction = probabilities[:,1]
    fpr, tpr, roc_threshold = roc_curve(test_y, prediction)
    auc = roc_auc_score(test_y, prediction)
    #print("Random Forest Classification AUC:", auc)
    return auc, fpr, tpr, roc_threshold, prediction

def getROC(expected, predicted):
    """
    Parameters
    ----------
    predictions : list of size len(complete dataset) of lists==[expected classifier label, prediction score of a given model, ...]

    Returns
    -------
    fpr : list
    tpr : list
    threshold : list
    auc : float
    """
    fpr, tpr, threshold = roc_curve(expected, predicted, pos_label=1)
    auc = roc_auc_score(expected, predicted)
    return fpr, tpr, threshold, auc

def getYoudensJ(fpr, tpr, threshold):
    """
    Parameters
    ----------
    fpr : list of false positive rates
    tpr : list of true positive rates
    threshold : threshold at calculated fpr, tpr

    Returns
    -------
    J : float, Youden's J statistic, where J = max(sensitivity + specificity -1)
    J_threshold : float, threshold of J --- commented out, but can be used if need be
    """
    J = 0.0
    #J_threshold = 0.0
    for t in range(len(threshold)):
        J_t = tpr[t] + (1-fpr[t]) - 1
        if J_t > J:
            J = J_t
            #J_threshold = threshold[t]
    #return J, J_threshold
    return J

def getMaxMCC(expected, predicted, threshold):
    """
    Parameters
    ----------
    expected : list of classifier labels
    predicted : list of a model's prediction values
    threshold : list of ROC thresholds to define boundary of P & N

    Returns
    -------
    maxMCC : float, maximum MCC of ROC curve
    """
    maxMCC = 0.0
    for t in enumerate(threshold):
        TP = 0
        FP = 0
        TN = 0
        FN = 0
        MCC = 0.0
        for i in range(len(expected)):
            if predicted[i] >= t[1] and expected[i] == 1:
                TP += 1
            if predicted[i] <= t[1] and expected[i] == 0:
                FP += 1
            if predicted[i] < t[1] and expected[i] == 1:
                FN += 1
            if predicted[i] < t[1] and expected[i] == 0:
                TN += 1
            numerator = (TP*TN) - (FP*FN)
            if (TP+FP)==0 or (TP+FN)==0 or (TN+FP)==0 or (TN+FN)==0:
                MCC = 0.0
            else:
                denominator = math.sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
                MCC = numerator / denominator
            if MCC > maxMCC:
                maxMCC = MCC
    return maxMCC

def summarizeStatistics(expected, predicted):
    """
    Parameters
    ----------
    expected : list of classifier labels
    predicted : list of a model's prediction values

    Returns
    -------
    fpr : list of false positive rate at given threshold
    tpr : list of true positive rate at given threshold
    threshold : list of arbitrary values where fpr or tpr changes 
    J : float, Youden's J statistic, where J = max(sensitivity + specificity -1)
    maxMCC : float, maximum Matthews Correlation Coefficient with given ROC
    """
    print("Getting summary statistics")
    fpr, tpr, threshold = roc_curve(expected, predicted, pos_label=1)
    auc = roc_auc_score(expected, predicted)
    J = getYoudensJ(fpr, tpr, threshold)
    maxMCC = getMaxMCC(expected, predicted, threshold)
    return fpr, tpr, threshold, auc, J, maxMCC

def compareROC(A_y, A_pred, B_y, B_pred):
    """
    Parameters
    ----------
    predA : list of lists [expected classifier label, prediction score of a given model, ...] of model A
    predB : list of lists [expected classifier label, prediction score of a given model, ...] of model B

    Returns
    -------
    p : float, Delong's p value
    """
    V_A10, V_A01 = delongs.structural_components(A_pred, A_y)
    V_B10, V_B01 = delongs.structural_components(B_pred, B_y)
    aucA = delongs.auc(A_pred, A_y)
    aucB = delongs.auc(B_pred, B_y)
    varA, varB, covarAB = delongs.computeCoVarMatrix(V_A10, V_A01, V_B10, V_B01, aucA, aucB) 
    z = delongs.z_score(varA, varB, covarAB, aucA, aucB)
    p = st.norm.sf(abs(z))*2
    return p

def main():

    # -------Set up input and output-------
    model = args.model
    label = args.label
    CBdist_byVirusPDB, R = getPDBdistances_byVirus(args.data, label) 
    features = []
    wts_byVirusPDB = [] 
    print("\nRunning AxIEM at upper boundary limit {} for {} epitope mapping predictions".format(R, label))
    wts_byVirusPDB = getNeighborWeights_byVirus(wts_byVirusPDB, CBdist_byVirusPDB, args.data, R, label ) 
    # Construct datasets:
    # 1) For each feature array, scale per-residue features {REU, CPrmsd, NV} using MinMaxScaler of range (0,1)
    # 2) Add Neighbor Sum feature using upper radius boundary value for scaled feature sets
    X, y, test_X, test_y = constructDatasets( args.data, CBdist_byVirusPDB, wts_byVirusPDB, label )
    AUC = 0.0
    fpr = []
    tpr = []
    thresh = []
    prd = []
    if model=='LR':
        AUC, fpr, tpr, thresh, prd = performLinearRegression(X, y, test_X, test_y)
        print("Linear regression AUC: ", AUC)
    if model=='BC':
        AUC, fpr, tpr, thresh, prd = performBayesClassification(X, y, test_X, test_y)
        print("Bayes classifier AUC: ", AUC)

    # Store new feature set values for plotting/analysis
    for i in range(len(test_X)):
        if model=='LR' or 'BC':
            features.append([str(test_y[i]), prd[i], str(test_X[i][0]), str(test_X[i][1]), str(test_X[i][2]), str(test_X[i][3])])
        if model=='BOTH':
            features.append([str(test_y[i]), lr_prd[i], bc_prd[i], str(test_X[i][0]), str(test_X[i][1]), str(test_X[i][2]), str(test_X[i][3])])

   # -------Write output for plotting--------
    print("Writing output")
    with open(args.feature_sets, 'w') as f:
        f.write("Classifier Prediction REU CPrmsd NV NS_R\n")
        for i in range(len(features)):
            for j in range(len(features[i])):
                f.write('%r ' % str(features[i][j]))
            f.write('\n')
            
    print("Done")

if __name__ == "__main__":
    main()
