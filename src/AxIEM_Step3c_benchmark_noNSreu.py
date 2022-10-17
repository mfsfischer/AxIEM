#!/usr/bin/env python3
"""
Performs leave-one-out analyses of epitope mapping models using a range of epitope radii
@author marionfischer
"""
# Parsing files
from argparse import ArgumentParser

# Basic calculations and statistical analysis
from collections import Counter, defaultdict
from math import sqrt
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
    "--randomized_data",
    dest="random",
    required=True,
    help="""INPUT randomized features file from Step 2 - see documentation for formatting"""
)
parser.add_argument(
    "--discotope",
    dest="discotope",
    required=True,
    help="""INPUT text file containing three columns [LABEL discotope_prediction VIRUS]"""
)
parser.add_argument(
    "--ellipro",
    dest="ellipro",
    required=True,
    help="""INPUT text file containing three columns [LABEL ellipro_prediction VIRUS]"""
)
parser.add_argument(
    "--summary",
    dest="summary",
    help="""OUTPUT file to write individual leave-one-out summary statistics"""
)
parser.add_argument(
    "--averages",
    dest="averages",
    help="""OUTPUT file to write mean and st. dev. of leave-one-out AUCS"""
)
parser.add_argument(
    "--rocs",
    dest="rocs",
    help="""OUTPUT file to write ROC curves for each leave-one-out test"""
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
    NSu_CP = pandas.DataFrame(data.iloc[:,u+1]).to_numpy()
    NSu_NV = pandas.DataFrame(data.iloc[:,u+2]).to_numpy()
    all_y = data.iloc[:,2].to_numpy()
    # Separate into training and test sets
    X = []
    y = []
    test_X = []
    test_y = []
    for v in range(len(viruses)):
        res_u_features = numpy.concatenate([perres_X[v],NSu_CP[v],NSu_NV[v]])
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

    print("\nStep 3 in AxIEM benchmark: Test model performance using leave-one-out tests")
    # For abbreviations: 'LR'=LinearRegression, 'BC'=BayesClassifier, 'LG'=LogisticRegression, 'RF'=RandomForestClassifier
    #                    'Random'=with randomized feature values

    # -------Set up------- 
    discotope_data = pandas.read_csv(args.discotope, delimiter=' ')
    disco_y = pandas.DataFrame(discotope_data.iloc[:,0]).to_numpy()
    disco_pred = pandas.DataFrame(discotope_data.iloc[:,1]).to_numpy()
    ellipro_data = pandas.read_csv(args.ellipro, delimiter=' ')
    elli_y = pandas.DataFrame(ellipro_data.iloc[:,0]).to_numpy()
    elli_pred = pandas.DataFrame(ellipro_data.iloc[:,1]).to_numpy()
    models = ['LR', 'BC', 'LG', 'RF', 'randLR', 'randBC', 'randLG', 'randRF']
    tests = getTestLabels(args.data)
    print(tests)
    # Initialize Output
    leave_out_AUC = {'8':{}, '16':{}, '24':{}, '32':{}, '40':{}, '48':{}, '56':{}, '64':{}, 'Discotope': {}, 'Ellipro': {}}
    roc_curves = []
    summarizeLR = {'8':{'mu':0.0,'sigma':0.0}, '16':{'mu':0.0,'sigma':0.0}, '24':{'mu':0.0,'sigma':0.0}, '32':{'mu':0.0,'sigma':0.0}, '40':{'mu':0.0,'sigma':0.0}, '48':{'mu':0.0,'sigma':0.0}, '56':{'mu':0.0,'sigma':0.0}, '64':{'mu':0.0,'sigma':0.0}}
    summarizeBC = {'8':{'mu':0.0,'sigma':0.0}, '16':{'mu':0.0,'sigma':0.0}, '24':{'mu':0.0,'sigma':0.0}, '32':{'mu':0.0,'sigma':0.0}, '40':{'mu':0.0,'sigma':0.0}, '48':{'mu':0.0,'sigma':0.0}, '56':{'mu':0.0,'sigma':0.0}, '64':{'mu':0.0,'sigma':0.0}}
    summarizeLG = {'8':{'mu':0.0,'sigma':0.0}, '16':{'mu':0.0,'sigma':0.0}, '24':{'mu':0.0,'sigma':0.0}, '32':{'mu':0.0,'sigma':0.0}, '40':{'mu':0.0,'sigma':0.0}, '48':{'mu':0.0,'sigma':0.0}, '56':{'mu':0.0,'sigma':0.0}, '64':{'mu':0.0,'sigma':0.0}}
    summarizeRF = {'8':{'mu':0.0,'sigma':0.0}, '16':{'mu':0.0,'sigma':0.0}, '24':{'mu':0.0,'sigma':0.0}, '32':{'mu':0.0,'sigma':0.0}, '40':{'mu':0.0,'sigma':0.0}, '48':{'mu':0.0,'sigma':0.0}, '56':{'mu':0.0,'sigma':0.0}, '64':{'mu':0.0,'sigma':0.0}}
    summarizeDE = {'Discotope':{'mu':0.0, 'sigma':0.0}, 'Ellipro':{'mu':0.0, 'sigma':0.0}}
    # -------Run each leave one out test by upper radius boundary and assess--------
    num_t = 0
    for t in tests:
        # Neighbor Sum values are dependent on upper boundary radius - NSu values for
        # each feature are stored in columns 3*i+3, u=upper radius and each NS feature
        # is ordered as 3+x, ordered as NSu_REU, NSu_CP, NSu_NV, and first three
        # columns are per-residue REU, CP, and NV features
        for i in range(8):
            upper_boundary = str((i+1)*8)
            X, y, test_X, test_y = constructDatasets(args.data, t, i*3+7)
            randX, y, test_randX, test_y = constructDatasets(args.random, t, i*3+7)
            # Create each leave-one-out model and assess
            print("With features")
            LR_AUC, LR_fpr, LR_tpr, LR_threshold, LR_prediction = performLinearRegression(X, y, test_X, test_y)
            BC_AUC, BC_fpr, BC_tpr, BC_threshold, BC_prediction = performBayesClassification(X, y, test_X, test_y)
            LG_AUC, LG_fpr, LG_tpr, LG_threshold, LG_prediction = performLogisticRegression(X, y, test_X, test_y)
            RF_AUC, RF_fpr, RF_tpr, RF_threshold, RF_prediction = performRandomForestClassification(X, y, test_X, test_y)
            print("AUCs:", LR_AUC, "\t", BC_AUC, "\t", LG_AUC,"\t", RF_AUC)
            # Store ROC curves for plotting
            for i in range(len(LR_fpr)):
                roc_curves.append([LR_fpr[i], LR_tpr[i], LR_threshold[i], t, upper_boundary, "LR"])
            for i in range(len(BC_fpr)):
                roc_curves.append([BC_fpr[i], BC_tpr[i], BC_threshold[i], t, upper_boundary, "BC"])
            for i in range(len(LG_fpr)):
                roc_curves.append([LG_fpr[i], LG_tpr[i], LG_threshold[i], t, upper_boundary, "LG"])
            for i in range(len(RF_fpr)):
                roc_curves.append([RF_fpr[i], RF_tpr[i], RF_threshold[i], t, upper_boundary, "RF"])
                
            # Perform negative control to ensure that feature set is providing a signal
            ## For each training set, randomize each feature's values based off of the
            ## mean and st. dev. protein's feature values - should expect an AUC of ~0.5
            print("With randomized features")
            randLR_AUC, randLR_fpr, randLR_tpr, randLR_threshold, randLR_prediction = performLinearRegression(randX, y, test_randX, test_y)
            randBC_AUC, randBC_fpr, randBC_tpr, randBC_threshold, randBC_prediction = performBayesClassification(randX, y, test_randX, test_y)
            randLG_AUC, randLG_fpr, randLG_tpr, randLG_threshold, randLG_prediction = performLinearRegression(randX, y, test_randX, test_y)
            randRF_AUC, randRF_fpr, randRF_tpr, randRF_threshold, randRF_prediction = performRandomForestClassification(randX,y,test_randX, test_y)
            print("AUCs:", randLR_AUC, "\t", randBC_AUC, "\t", randLG_AUC,"\t", randRF_AUC)

            # Individual leave-one-out AUC statistics by method using upper boundary and by virus
            leave_out_AUC[upper_boundary][t] = {}
            leave_out_AUC[upper_boundary][t]['LR']      = LR_AUC
            leave_out_AUC[upper_boundary][t]['BC']      = BC_AUC
            leave_out_AUC[upper_boundary][t]['LG']      = LG_AUC
            leave_out_AUC[upper_boundary][t]['RF']      = RF_AUC
            leave_out_AUC[upper_boundary][t]['randLR']  = randLR_AUC
            leave_out_AUC[upper_boundary][t]['randBC']  = randBC_AUC
            leave_out_AUC[upper_boundary][t]['randLG']  = randLG_AUC
            leave_out_AUC[upper_boundary][t]['randRF']  = randRF_AUC

            # Mean and st. dev. of upper boundary AUC performance in one pass
            summarizeLR[upper_boundary]['mu'] += LR_AUC
            summarizeLR[upper_boundary]['sigma'] += LR_AUC * LR_AUC
            summarizeBC[upper_boundary]['mu'] += BC_AUC
            summarizeBC[upper_boundary]['sigma'] += BC_AUC * BC_AUC
            summarizeLG[upper_boundary]['mu'] += LG_AUC
            summarizeLG[upper_boundary]['sigma'] += LG_AUC * LG_AUC
            summarizeRF[upper_boundary]['mu'] += RF_AUC
            summarizeRF[upper_boundary]['sigma'] += RF_AUC * RF_AUC
            print("Finished running leave out {} tests for boundary {}\n".format(t, upper_boundary))
            num_t += 1
    print("\nBenchmarking Discotope and Ellipro performance")
    # Individual virus prediction performance of Discotope
    discotope_AUC = 0.0
    disco_y = []
    disco_pred = []
    for t in tests:
        if discotope_data['Virus'].eq(t.strip()).any():
            t_discotope = discotope_data[discotope_data.Virus == t]
            disco_y  = t_discotope['Classifier'].values
            disco_pred = t_discotope['Prediction'].values
            t_fpr, t_tpr, t_threshold, t_AUC = getROC( disco_y, disco_pred )
            leave_out_AUC['Discotope'][t] = t_AUC
            summarizeDE['Discotope']['mu'] += t_AUC
            summarizeDE['Discotope']['sigma'] += t_AUC * t_AUC
            for i in range(len(t_fpr)):
                roc_curves.append([t_fpr[i], t_tpr[i], t_threshold[i], t, "Discotope", "Discotope"])
        else:
            print(t, "has no Discotope predictions")
    # Individual virus prediction performance of Ellipro
    ellipro_AUC = 0.0
    elli_y = []
    elli_pred = []
    for t in tests:
        if ellipro_data['Virus'].eq(t).any():
            t_ellipro = ellipro_data[ellipro_data.Virus == t]
            elli_y  = t_ellipro['Classifier'].values
            elli_pred = t_ellipro['Prediction'].values
            t_fpr, t_tpr, t_threshold, t_AUC = getROC( elli_y, elli_pred )
            leave_out_AUC['Ellipro'][t] = t_AUC
            summarizeDE['Ellipro']['mu'] += t_AUC
            summarizeDE['Ellipro']['sigma'] += t_AUC * t_AUC
            for i in range(len(t_fpr)):
                roc_curves.append([t_fpr[i], t_tpr[i], t_threshold[i], t, "Ellipro", "Ellipro"])
        else:
            print(t, "has no Ellipro predictions")

    print("Writing output")
    # -------Fig2A--------
    with open(args.summary, 'w') as ind:
        ind.write("Boundary Virus Model AUC\n")
        for u in ['8','16','24','32','40','48','56','64']:
            for t in tests:
                for m in models:
                    ind.write('%r %r %r %r\n' % (u, t, m, leave_out_AUC[u][t][m]))
        for d in ['EBOV','H3','H7','HIV','RSV']:
            ind.write('Discotope %r Discotope %r\n' % (d, leave_out_AUC['Discotope'][d]))
        for t in tests:
            ind.write('Ellipro %r Ellipro %r\n' % (t, leave_out_AUC['Ellipro'][t]))
    # -------Fig2B: Finalize Mean and St. Dev. calculations-------
    for u in ['8','16','24','32','40','48','56','64']:
        summarizeLR[u]['mu'] /= num_t
        summarizeBC[u]['mu'] /= num_t
        summarizeLG[u]['mu'] /= num_t
        summarizeRF[u]['mu'] /= num_t
        summarizeLR[u]['sigma'] = sqrt(summarizeLR[u]['sigma'] / (num_t-summarizeLR[u]['mu']))
        summarizeBC[u]['sigma'] = sqrt(summarizeBC[u]['sigma'] / (num_t-summarizeBC[u]['mu']))
        summarizeLG[u]['sigma'] = sqrt(summarizeLG[u]['sigma'] / (num_t-summarizeLG[u]['mu']))
        summarizeRF[u]['sigma'] = sqrt(summarizeRF[u]['sigma'] / (num_t-summarizeRF[u]['mu']))
    summarizeDE['Discotope']['mu'] /= num_t
    summarizeDE['Ellipro']['mu'] /= num_t
    summarizeDE['Discotope']['sigma'] = sqrt(summarizeDE['Discotope']['sigma'] / (num_t-summarizeDE['Discotope']['mu'])) 
    summarizeDE['Ellipro']['sigma'] = sqrt(summarizeDE['Ellipro']['sigma'] / (num_t-summarizeDE['Ellipro']['mu']) )
    with open(args.averages, 'w') as avg:
        avg.write("Boundary Model MeanAUC stDevAUC\n")
        for u in ['8','16','24','32','40','48','56','64']:
            avg.write('%r LR %r %r\n' % (u, summarizeLR[u]['mu'], summarizeLR[u]['sigma']))
            avg.write('%r BC %r %r\n' % (u, summarizeBC[u]['mu'], summarizeBC[u]['sigma']))
            avg.write('%r LG %r %r\n' % (u, summarizeLG[u]['mu'], summarizeLG[u]['sigma']))
            avg.write('%r RF %r %r\n' % (u, summarizeRF[u]['mu'], summarizeRF[u]['sigma']))
        avg.write('Discotope Discotope %r %r\n' % (summarizeDE['Discotope']['mu'], summarizeDE['Discotope']['sigma']))
        avg.write('Ellipro Ellipro %r %r\n' % (summarizeDE['Ellipro']['mu'], summarizeDE['Ellipro']['sigma']))
    # -------Supp Info: ROC curves for each leave out test-------
    with open(args.rocs, 'w') as roc:
        roc.write("tpr fpr threshold Virus Boundary Model\n")
        for i in range(len(roc_curves)):
            roc.write(" ".join(map(str, roc_curves[i])))
            roc.write("\n")
    print("Done")

if __name__ == "__main__":
    main()
