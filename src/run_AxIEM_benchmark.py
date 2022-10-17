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
    """Benchmarks AxIEM appproach by:
    1) Runs leave-one-out predictions for a range of indicated radii and generates summary statistics for each leave-one-out test,
    2) Calculates statistical analysis for all test cases with the same epitope upper boundary radius used to calculate the Neighbor Sum feature,
    3) Compares performance of model type (e.g. Bayes classifier or random forest classifier) to Discotope and Ellipro predictions"""
    )
parser.add_argument(
    "--data",
    dest="data",
    required=True,
    help="""INPUT text file complete_dataset.txt"""
    )
parser.add_argument(
    "--discotope",
    dest="discotope",
    required=True,
    help="""INPUT text file containing two columns [classifier_label discotope_prediction]"""
    )
parser.add_argument(
    "--ellipro",
    dest="ellipro",
    required=True,
    help="""INPUT text file containing two columns [classifier_label ellipro_prediction]"""
    )
parser.add_argument(
    "--individual_summary",
    dest="individual_summary",
    help="""OUTPUT file to write individual leave-one-out summary statistics"""
)
parser.add_argument(
    "--features",
    dest="feature_sets",
    help="""OUTPUT file to write feature sets including raw Neighbor Sum values and scaled feature sets"""
)

# Parse the command line
args = parser.parse_args()
# Parse PDB file
pdbparser = PDBParser()

def getPDBdistances_byVirus(dataset):
    """
    Parameters
    ----------
    dataset : txt file of annotated data (e.g. AxIEM.data)

    Returns
    -------
    viruses_cb_distances : dict of dict containing numpy arrays size (n,n) of each virus's PDB's contact map distances 
    """
    data = pandas.read_csv(dataset, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    pdbs = data.iloc[:,1].to_list()
    viruses_cb_distances = defaultdict(lambda : defaultdict([]))
    for i in range(len(viruses)):
        if not viruses[i] in viruses_cb_distances:
            viruses_cb_distances[viruses[i]] = {}
            if not pdbs[i] in viruses_cb_distances[viruses[i]]:
                pdb = "pdb_structures/{}".format(pdbs[i])
                structure = pdbparser.get_structure("apdb", pdb)
                residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
                print("Getting {} structure {} with {} residues".format(viruses[i], pdbs[i], len(residues)))
                viruses_cb_distances[viruses[i]][pdbs[i]] = getCbDistances(residues)
        else:
            if not pdbs[i] in viruses_cb_distances[viruses[i]]:
                pdb = "pdb_structures/{}".format(pdbs[i])
                structure = pdbparser.get_structure("apdb", pdb)
                residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
                print("Getting {} structure {} with {} residues".format(viruses[i], pdbs[i], len(residues)))
                viruses_cb_distances[viruses[i]][pdbs[i]] = getCbDistances(residues)
    return viruses_cb_distances

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

def getNeighborWeights_byVirus(wts, cb_distances, dataset, upper ):
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
    X                : numpy array size (len(n-m), 4) of raw training feature set values + NS_u calculated with raw data
    Xscaled          : numpy array list size (len(n-m), 4) of scaled training feature set + NS_u calculated with scaled data
    randX            : numpy array size (len(n-m), 4) of randomized raw training feature set values + NS_u calculated with randomized raw data 
    randXscaled      : numpy array list size (len(n-m), 4) of randomized scaled training feature set + NS_u calculated with randomized scaled data
    y                : numpy array size len(m) of classifier labels
    test_X           : numpy array size (len(m), 4) of raw test feature set values + NS_u calculated with raw data
    test_Xscaled     : numpy array list size (len(m), 4) of scaled test feature set + NS_u calculated with scaled data
    test_randX       : numpy array size (len(m), 4) of randomized raw test feature set values + NS_u calculated with randomized raw data 
    test_randXscaled : numpy array list size (len(m), 4) of randomized scaled test feature set + NS_u calculated with randomized scaled data
    test_y           : numpy array size len(m) of classifier labels
    """
    data = pandas.read_csv(data, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    pdbs = data.iloc[:,1].to_list()
    allX = pandas.DataFrame(data.iloc[:,4:7]).to_numpy()
    scaler = MinMaxScaler()
    allXscaled = scaler.fit_transform(allX)
    allY = data.iloc[:,2].to_numpy()

    X = []
    Xscaled = []
    randX = []
    randXscaled = []
    y = []

    test_X = []
    test_Xscaled = []
    test_randX = []
    test_randXscaled = []
    test_y = []

    # Add Neighbor Sum (NS_u) feature using raw or scaled per-residue features
    allX, allXscaled = assignNeighborSum( allX, allXscaled, pdbs, viruses, virus_cb_distances, virus_wts)
    # Randomize features for negative control testing
    allX_rand = randomize(allX)
    allXscaled_rand = randomize(allXscaled)
    # Separate into training and test sets
    for v in range(len(viruses)):
        if viruses[v] == test_protein:
            test_X.append(allX[v])
            test_Xscaled.append(allXscaled[v])
            test_randX.append(allX_rand[v])
            test_randXscaled.append(allXscaled[v])
            test_y.append(allY[v])
        else:
            X.append(allX[v])
            Xscaled.append(allXscaled[v])
            randX.append(allX_rand[v])
            randXscaled.append(allXscaled_rand[v])
            y.append(allY[v])
    return numpy.asarray(X), numpy.asarray(Xscaled), numpy.asarray(randX), numpy.asarray(randXscaled), numpy.asarray(y), numpy.asarray(test_X), numpy.asarray(test_Xscaled), numpy.asarray(test_randX), numpy.asarray(test_randXscaled), numpy.asarray(test_y)

def assignNeighborSum( X_in, X_in_scaled, pdbs, viruses, CBdicts, WTSdicts ):
    """
    Parameters
    ----------
    X_in : numpy array size (n, in) where in==[REU, CPrmsd, NV]
    X_in_scaled : numpy array size (n, scaled) where scaled==MinMaxScaler([REU, CPrmsd, NV], axis=0)
    pdbs : list size (n) of pdb names
    viruses : list size (n) of virus names
    CBdicts : dict of dicts to get numpy array size (n,n) of cb_distances for each PDB
    upper : float value to use as upper boundary for NS_u feature calculated in getNeighborWeights

    Returns
    -------
    X : numpy array size (n, p) where p==[REU, CPrmsd, NV, NS_u]
    X_scaled : numpy array size (n, p) where p==[scaled-REU, scaled-CPrmsd, scaled-NV, NS_u with scaled features], axis=0)

    """
    # -----Extract residue predictors by PDB------
    pdb_counts = dict(Counter(pdbs).items())         # Get number, lengths, names of PDBs in dataset - all ordered
    num_pdbs = len(pdb_counts)                       # Number of PDBs
    pdb_lengths = list(pdb_counts.values())          # Lengths of PDBs
    pdb_names = list(pdb_counts.keys())              # Names of PDBs
    
    pdb_index = 0
    res_index = 0
    X = numpy.zeros((X_in.shape[0], X_in.shape[1]+1))
    X_scaled = numpy.zeros((X_in_scaled.shape[0], X_in_scaled.shape[1]+1))
    # For each PDB in the dataset
    for p in range(num_pdbs):
        pdb_data = []
        pdb_data_scaled = []
        pdb = "pdb_structures/{}".format(pdb_names[p])
        virus = viruses[res_index]
        pdb_cb_dists = CBdicts[virus][pdb_names[p]]
        pdb_wts = WTSdicts[virus][pdb_names[p]]
        # Get all residues within single PDB 
        for res in range(pdb_lengths[p]):
            pdb_data.append(X_in[pdb_index])
            pdb_data_scaled.append(X_in_scaled[pdb_index])
            pdb_index += 1
        # Add Neighbor Sum feature given weights calculated with a given upper boundary
        pdb_NS, pdb_NS_scaled = addNeighborSum(pdb_cb_dists, pdb_wts, numpy.asarray(pdb_data), numpy.asarray(pdb_data_scaled))
        # Reconstruct entire data set with same indexing as X_in/X_in_scaled
        for res in range(pdb_NS.shape[0]):
            X[res_index] = pdb_NS[res]
            X_scaled[res_index] = pdb_NS_scaled[res]
            res_index += 1
    return X, X_scaled

def addNeighborSum(cb_distances, wts, features, scaled_features):
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
    # Find inner sum of all per-residue features/scaled features
    inner_sum = numpy.zeros(features.shape[0])
    inner_sum_scaled = numpy.zeros(features.shape[0])
    for res in range(features.shape[0]):
        for p in range(features.shape[1]):
            inner_sum[res] += features[res][p]
            inner_sum_scaled[res] += scaled_features[res][p]
    # Find outer sum of all surrounding weighted residues' features
    # accounting for all distances to each residue
    # and return new features as 
    # {f1, f2, ..., fn, neighbor_score} or
    # {f1, f2, ..., fn, neighbor_scaled_score}
    X = numpy.zeros( (features.shape[0], features.shape[1]+1) )
    X_scaled = numpy.zeros( (features.shape[0], features.shape[1]+1) )
    for res in range(len(inner_sum)):
        # Get weighted neighborhood weights for reference residue
        neighbor_score = 0
        neighbor_scaled_score = 0
        for w in range(len(wts)):
            neighbor_score += wts[res][w] * inner_sum[res]
            neighbor_scaled_score += wts[res][w] * inner_sum_scaled[res]
        for p in range(features.shape[1]+1):
            if p < features.shape[1]:
                X[res][p] = features[res][p]
                X_scaled[res][p] = scaled_features[res][p]
            if p == features.shape[1]:
                X[res][p] = neighbor_score
                X_scaled[res][p] = neighbor_scaled_score
    return X, X_scaled

def randomize( features ):
    """
    Parameters
    ----------
    features : numpy array size (n, 4) of input feature value array

    Returns
    -------
    randomized_features = numpy array size (n, 4) of features that have been randomized using the gaussian distribution's mu and sigma of each feature
    """
    randomized_features = numpy.zeros((features.shape[0], features.shape[1]))
    mu = numpy.mean(features, axis=0)
    sigma = numpy.std(features, axis=0)
    for residue in range(features.shape[0]):
        for feature in range(features.shape[1]):
            randomized_features[residue][feature] = random.gauss(mu[feature], sigma[feature])
    return randomized_features

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
    print("Starting benchmark")
    # For abbreviations, 'LR'=LinearRegression, 'BC'=BayesClassifier, 'RF'=RandomForestClassifier,
    # 'S-' or 's'= MinMaxScaled per-residue feature values, 'Random'=randomized feature values
    #models = ['LinearRegression', 'BayesClassifier', 'RandomForest', 'RandomLR', 'RandomBC', 'RandomRF','S-LinearRegression', 'S-BayesClassifier', 'S-RandomForest', 'S-RandomLR', 'S-RandomBC', 'S-RandomRF', 'Discotope', 'Ellipro']
    models = ['LR', 'BC', 'RF', 'sLR', 'sBC', 'sRF', 'randLR', 'randBC', 'randRF', 'srandLR', 'srandBC', 'srandRF']
    proteins = ['EBOV','H3','H7','HIV','RSV','SARS','SARS2']
    prot = ['EBOV', 'H3', 'H7']
    protein_len = [708, 339, 1422, 1339, 1236, 2970, 2706]
    prot_len = [708, 339]
    radii = ['8','16','24','32','40','48','56','64','R']
    rad = ['32','R']
    abv_models = ['LR', 'BC', 'RF', 'sLR', 'sBC', 'sRF']
    discotope_data = pandas.read_csv(args.discotope, delimiter=' ')
    disco_y = pandas.DataFrame(discotope_data.iloc[:,0]).to_numpy()
    disco_pred = pandas.DataFrame(discotope_data.iloc[:,1]).to_numpy()
    ellipro_data = pandas.read_csv(args.ellipro, delimiter=' ')
    elli_y = pandas.DataFrame(ellipro_data.iloc[:,0]).to_numpy()
    elli_pred = pandas.DataFrame(ellipro_data.iloc[:,1]).to_numpy()
    CBdist_byVirusPDB = getPDBdistances_byVirus(args.data) 
    # Output
    leave_out_AUC = {'8':{}, '16':{}, '24':{}, '32':{}, '40':{}, '48':{}, '56':{}, '64':{}, 'R':{}, 'Discotope': {}, 'Ellipro': {}}
    summarize = {'8':{}, '16':{}, '24':{}, '32':{}, '40':{}, '48':{}, '56':{}, '64':{}, 'R':{}, 'Discotope': {'AUC':0.0, 'J':0.0, 'maxMCC':0.0}, 'Ellipro': {'AUC':0.0, 'J':0.0, 'maxMCC':0.0}}
    #summary_statistics = ['AUC', 'J', 'maxMCC']
    #summarize = dict.fromkeys(abv_models, dict.fromkeys(radii, dict.fromkeys(summary_statistics, 0.0)))
    model_ROC = []

    # -------Run each leave one out test by upper radius boundary and assess--------
    # Get feature sets with upper boundary-dependent Neighbor Sum feature + scaled feature values
    features = []
    scaled_features = []
    for u in radii:
        wts_byVirusPDB = [] 
        if u != 'R':
            wts_byVirusPDB = getNeighborWeights_byVirus(wts_byVirusPDB, CBdist_byVirusPDB, args.data, float(u) )
        pro = 0
        for p in proteins:
            if u != 'R':
                print("\nPerforming leave-one-out tests using {} at upper boundary limit {}".format(p, u))
            # When the test case's protein length is used to determine the upper boundary
            else:
                upper_radius = 15.1 + 0.0147 * protein_len[pro]
                print("\nPerforming leave-one-out tests using {} at upper boundary limit {}".format(p, upper_radius))
                wts_byVirusPDB = getNeighborWeights_byVirus(wts_byVirusPDB, CBdist_byVirusPDB, args.data, upper_radius ) 
            # Construct datasets:
            # 1) For each feature array, scale per-residue features {REU, CPrmsd, NV} using MinMaxScaler of range (0,1)
            # 2) Add Neighbor Sum feature using upper radius boundary value for raw and scaled feature sets
            # 3) Randomize feature set values by gaussian(mu, sigma) of feature's values
            X, Xscaled, randX, randXscaled, y, test_X, test_Xscaled, test_randX, test_randXscaled, test_y = constructDatasets( args.data, CBdist_byVirusPDB, wts_byVirusPDB, p )

            # Create each leave-one-out model and assess
            print("With features - not scaled")
            LR_AUC, LR_fpr, LR_tpr, LR_threshold, LR_prediction = performLinearRegression(X, y, test_X, test_y)
            BC_AUC, BC_fpr, BC_tpr, BC_threshold, BC_prediction = performBayesClassification(X, y, test_X, test_y)
            RF_AUC, RF_fpr, RF_tpr, RF_threshold, RF_prediction = performRandomForestClassification(X, y, test_X, test_y)
            print("AUC:", LR_AUC, "\t", BC_AUC, "\t", RF_AUC)
            
            print("With features - scaled")
            sLR_AUC, sLR_fpr, sLR_tpr, sLR_threshold, sLR_prediction = performLinearRegression(Xscaled, y, test_Xscaled, test_y)
            sBC_AUC, sBC_fpr, sBC_tpr, sBC_threshold, sBC_prediction = performBayesClassification(Xscaled, y, test_Xscaled, test_y)
            sRF_AUC, sRF_fpr, sRF_tpr, sRF_threshold, sRF_prediction = performRandomForestClassification(Xscaled, y, test_Xscaled, test_y)
            print("AUCs:", sLR_AUC, sBC_AUC, sRF_AUC)
            
            # Perform negative control to ensure that feature set is providing a signal
            ## For each training set, randomize each feature's values based off of the
            ## mean and st. dev. protein's feature values - should expect an AUC of ~0.5
            print("With randomized features - not scaled")
            randLR_AUC, randLR_fpr, randLR_tpr, randLR_threshold, randLR_prediction = performLinearRegression(randX, y, test_randX, test_y)
            randBC_AUC, randBC_fpr, randBC_tpr, randBC_threshold, randBC_prediction = performBayesClassification(randX, y, test_randX, test_y)
            randRF_AUC, randRF_fpr, randRF_tpr, randRF_threshold, randRF_prediction = performRandomForestClassification(randX,y,test_randX, test_y)
            print("AUCs:", randLR_AUC, randBC_AUC, randRF_AUC)
            
            print("With randomized features - scaled")
            srandLR_AUC, srandLR_fpr, srandLR_tpr, srandLR_threshold, srandLR_prediction = performLinearRegression(randXscaled, y, test_randXscaled, test_y)
            srandBC_AUC, srandBC_fpr, srandBC_tpr, srandBC_threshold, srandBC_prediction = performBayesClassification(randXscaled, y, test_randXscaled, test_y)
            srandRF_AUC, srandRF_fpr, srandRF_tpr, srandRF_threshold, srandRF_prediction = performRandomForestClassification(randXscaled, y, test_randXscaled, test_y)
            print("AUCs:", srandLR_AUC, srandBC_AUC, srandRF_AUC)
            # Individual leave-one-out AUC statistics by method using upper boundary and by virus
            leave_out_AUC[u][p] = {}
            leave_out_AUC[u][p]['LR']      = LR_AUC
            leave_out_AUC[u][p]['BC']      = BC_AUC
            leave_out_AUC[u][p]['RF']      = RF_AUC
            leave_out_AUC[u][p]['randLR']  = randLR_AUC
            leave_out_AUC[u][p]['randBC']  = randBC_AUC
            leave_out_AUC[u][p]['randRF']  = randRF_AUC
            leave_out_AUC[u][p]['sLR']     = sLR_AUC
            leave_out_AUC[u][p]['sBC']     = sBC_AUC
            leave_out_AUC[u][p]['sRF']     = sRF_AUC
            leave_out_AUC[u][p]['srandLR'] = srandLR_AUC
            leave_out_AUC[u][p]['srandBC'] = srandBC_AUC
            leave_out_AUC[u][p]['srandRF'] = srandRF_AUC
            # Store new feature set values for plotting/analysis
            for i in range(len(test_X)):
                features.append([str(test_y[i]), str(test_X[i][0]), str(test_X[i][1]), str(test_X[i][2]), str(test_X[i][3]), p, u])
                scaled_features.append([str(test_y[i]), str(test_Xscaled[i][0]), str(test_Xscaled[i][1]), str(test_Xscaled[i][2]), str(test_Xscaled[i][3]), p, u])
            pro += 1
        print("Finished running leave out bests for boundary {}".format(u))

    print("\nBenchmarking Discotope and Ellipro performance")
    # Individual virus prediction performance of Discotope
    discotope_AUC = 0.0
    disco_y = []
    disco_pred = []
    for p in proteins:
        if discotope_data['Virus'].eq(p).any():
            p_discotope = discotope_data[discotope_data.Virus == p]
            disco_y  = p_discotope['Classifier'].values
            disco_pred = p_discotope['Prediction'].values
            p_fpr, p_tpr, p_threshold, p_AUC = getROC( disco_y, disco_pred )
            leave_out_AUC['Discotope'][p] = p_AUC
        else:
            print(p, "has no Discotope predictions")
    # Individual virus prediction performance of Ellipro
    ellipro_AUC = 0.0
    elli_y = []
    elli_pred = []
    for p in proteins:
        if ellipro_data['Virus'].eq(p).any():
            p_ellipro = ellipro_data[ellipro_data.Virus == p]
            elli_y  = p_ellipro['Classifier'].values
            elli_pred = p_ellipro['Prediction'].values
            p_fpr, p_tpr, p_threshold, p_AUC = getROC( elli_y, elli_pred )
            leave_out_AUC['Ellipro'][p] = p_AUC
        else:
            print(p, "has no Ellipro predictions")

   # -------Write output for plotting--------
    print("Writing output")
    with open(args.individual_summary, 'w') as ind, open(args.feature_sets, 'w') as f:
        # Data for Fig2A
        # protein_AUC: [radii][proteins][models]
        ind.write("Boundary Virus Model AUC\n")
        for u in radii:
            for p in proteins:
                for m in models:
                    ind.write('%r %r %r %r\n' % (u, p, m, leave_out_AUC[u][p][m]))
        for d in ['EBOV','H3','H7','HIV','RSV']:
            ind.write('Discotope %r Discotope %r\n' % (d, leave_out_AUC['Discotope'][d]))
        for p in proteins:
            ind.write('Ellipro %r Ellipro %r\n' % (p, leave_out_AUC['Ellipro'][p]))

        # Data for feature distributions (Supplementary Fig1+Fig2, overlap n scores)
        f.write("Classifier REU CPrmsd NV NS_u Protein u_Boundary Scaler\n")
        for i in range(len(features)):
            for j in range(len(features[i])):
                f.write('%r ' % str(features[i][j]))
            f.write('None \n')
        for i in range(len(scaled_features)):
            for j in range(len(scaled_features[i])):
                f.write('%r ' % str(scaled_features[i][j]))
            f.write('MinMax \n')
            
    print("Done")

if __name__ == "__main__":
    main()
