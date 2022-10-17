
# AxIEM

## Description
Antibody epitope mapping of viral proteins plays a vital role in understanding immune system mechanisms of protection. In the case of class I viral fusion proteins, recent advances in cryo-electron microscopy and protein stabilization techniques have highlighted the importance of cryptic or ‘alternative’ conformations that expose epitopes targeted by potent neutralizing antibodies. Thorough epitope mapping of such metastable conformations is difficult but is critical for understanding sites of vulnerability in class I fusion proteins that occur as transient conformational states during viral attachment and fusion. We introduce a novel method Accelerated class I fusion protein Epitope Mapping (AxIEM) that accounts for fusion protein flexibility to improve out-of-sample prediction of discontinuous antibody epitopes. Harnessing data from previous experimental epitope mapping efforts of several class I fusion proteins, we demonstrate that accuracy of epitope prediction depends on residue environment and allows for the prediction of conformation-dependent antibody target residues. We also show that AxIEM can identify common epitopes and provide structural insights for the development and rational design of vaccines.

email: marion.f.s.fischer@gmail.com

## Description of benchmark

Please refer to the [protocol capture](Protocol_Capture.pdf) for details. In brief the benchmark was run as follows.

The initial annotated dataset `AxIEM.data` was constructed by pasting and concatenating the classifier labels and pre-computed Rosetta per-residue total score energies for all PDB structures, to which virus protein name, PDB ID, and PDB residue IDs labels were added for data clarity. Next, the contact proximity variation and neighbor vector features were calculated using the following script.
```
python src/AxIEM_Step1_benchmark.py --data AxIEM.data --features AxIEM_per-residue.features
```

Afterwards, Neighbor Sums were calculated and appended to generate the complete dataset `AxIEM_updated.features`. Randomized feature values can be found in `AxIEM_randomized.features`.
```
python src/AxIEM_Step2_benchmark.py --data AxIEM_per-residue.features
--features AxIEM_updated.features --randomized_features AxIEM_randomized.features \
--plotting feature_distributions.txt
```

Finally, linear regression, Bayes classifier, Logistic regression, and random forest classifier models were trained and test using leave-out tests. Individual leave-out performance tests can be found in the `results/benchmark_all_leaveout_AUCs`.
```
python src/AxIEM_Step3_benchmark.py --data AxIEM_updated.features \
--randomized_data AxIEM_randomized.features --discotope Discotope.data \
--ellipro Ellipro.data --summary benchmark_leaveout_AUCs.txt \
--averages benchmark_avgAUC.txt --rocs benchmark_rocs.txt
```
