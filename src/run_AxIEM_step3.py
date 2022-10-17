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
import pomegranate
from pomegranate import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

# Create command line parser
parser = ArgumentParser(
    """Benchmarks AxIEM appproach by:
    1) Calculates per-residue features CP (contact proximity) and NV (Neighbor Vector) using at least two PDB structures
       of identical protein length. PDBs should contain similar/same virus protein sequence.
    2) Computes 'neighborhood' feature scores using indicated upper boundary (u)
    3) Runs leave-one-out predictions for each protein and generates summary statistics for each leave-one-out test
    """
)
parser.add_argument(
    "--data",
    dest="data",
    required=True,
    help="""INPUT features file from Step 2 - see documentation for formatting"""
)
parser.add_argument(
    "--summary",
    dest="summary",
    help="""OUTPUT file to write individual leave-one-out summary statistics"""
)
parser.add_argument(
    "--virus",
    dest="virus"
)

# Parse the command line
args = parser.parse_args()

def getTestLabels(dataset):
    """
    Parameters
    ----------
    dataset :

    Returns
    -------
    test_labels :
    """
    data = pandas.read_csv(dataset, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    test_labels = list(numpy.unique(viruses))
    test_labels.pop(0)
    return test_labels

def constructDatasets(dataset, test_protein, u):
    """
    Parameters
    ----------
    dataset      : text file of annotated complete dataset with each line [VIRUS PDB CLASSIFIER_LABEL RES_ID REU CPrmsd NV] of len(n)
    test_protein : string of protein name to be left out of training set and to be used for testing
    u            : int, index of NS_u REU in dataset

    Returns
    -------
    X            : numpy array size (len(n-m), 30) of training feature set values
    randX        : numpy array size (len(n-m), 30) of randomized training feature set values
    y            : numpy array size len(n-m) of classifier labels
    test_X       : numpy array size (len(m), 30) of test feature set values
    test_randX   : numpy array size (len(m), 30) of randomized test feature set values
    test_y       : numpy array size len(m) of classifier labels
    """
    data = pandas.read_csv(dataset, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    pdbs = data.iloc[:,1].to_list()
    perres_X = pandas.DataFrame(data.iloc[:,4:7]).to_numpy()
    NSu_X = pandas.DataFrame(data.iloc[:,u:u+3]).to_numpy()
    all_y = data.iloc[:,2].to_numpy()
    # Separate into training and test sets
    X = []
    y = []
    test_X = []
    test_y = []
    for v in range(len(viruses)):
        res_u_features = numpy.concatenate([perres_X[v],NSu_X[v]])
        if viruses[v] == test_protein:
            test_X.append(res_u_features)
            test_y.append(all_y[v])
        else:
            # Exclude fusion proteins if related to test set's protein:
            if test_protein=='H3' and viruses[v]!='H7':
                X.append(res_u_features)
                y.append(all_y[v])
            elif test_protein=='H7' and viruses[v]!='H3':
                X.append(res_u_features)
                y.append(all_y[v])
            if test_protein=='SARS' and viruses[v]!='SARS2':
                X.append(res_u_features)
                y.append(all_y[v])
            if test_protein=='SARS2' and viruses[v]!='SARS':
                X.append(res_u_features)
                y.append(all_y[v])
            else:
                X.append(res_u_features)
                y.append(all_y[v])
    return numpy.asarray(X), numpy.asarray(y), numpy.asarray(test_X), numpy.asarray(test_y)

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

def performLogisticRegression(train_X, train_y, test_X, test_y):
    model = LogisticRegression().fit(train_X, train_y)
    prediction = model.predict(test_X)
    fpr, tpr, roc_threshold = roc_curve(test_y, prediction)
    auc = roc_auc_score(test_y, prediction)
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

def getAccuracyLabels(expected, predicted, maxMCC):
    labels = numpy.empty(len(predicted))
    print(maxMCC)
    for i in range(len(predicted)):
        if predicted[i] >= maxMCC and expected[i] == 1:
            print(predicted[i])
            labels[i] = 0
        if predicted[i] <= maxMCC and expected[i] == 0:
            labels[i] = 1
        if predicted[i] < maxMCC and expected[i] == 1:
            labels[i] = 2
        if predicted[i] < maxMCC and expected[i] == 0:
            labels[i] = 3
    return labels

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

def main():

    # For abbreviations: 'LR'=LinearRegression, 'BC'=BayesClassifier, 'LG'=LogisticRegression, 'RF'=RandomForestClassifier
    #                    'Random'=with randomized feature values

    virus = args.virus
    # -------Build AxIEM and assess--------
    X, y, test_X, test_y = constructDatasets(args.data, args.virus, 3*3+7)
    LR_AUC, LR_fpr, LR_tpr, LR_threshold, LR_prediction = performLinearRegression(X, y, test_X, test_y)
    print("AUC:", LR_AUC)

    # -------Assess for TP, TN, FP, FN---------
    J = getYoudensJ(test_y, LR_prediction, LR_threshold)
    labels = getAccuracyLabels(test_y, LR_prediction, J)
    with open(args.summary, 'w') as lab:
        for i in range(len(labels)):
            lab.write(str(labels[i])+"\n")
    print("Done")

if __name__ == "__main__":
    main()
