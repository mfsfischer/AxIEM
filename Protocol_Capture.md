---
title: Protocol capture
---

# Preparation of Input PDB files for Feature Calculations

## Initial Curation of aligned PDB files

Class I fusion proteins are trimeric, which may or may not share C3
symmetry, and undergo multiple Ångstrom structural rearrangements during
cellular entry. The features used to build the AxIEM model require that
all conformations, *i.e.* For all methods described below, all PDB
structures used to represent a class I fusion protein during cellular
entry must share the same chain and residue order with respect to the
superimposition of all structures of the same fusion protein. To meet
this criteria, each viral fusion protein used in this study must have at
least two PDB structures within the [Protein Data Bank](www.rcsb.org)
that share no less than 1.0 Å RMSD to any other determined structure,
and that the PDBs used for the dataset share congruous chain and
sequence order with respect to the superimposed PDB ensemble.\
Curation of structures began with a PDB search of known class I fusion
proteins in their trimeric state, narrowed down by *i)* virus family,
*ii)* strain or serotype with at least two structures, and *iii)* all
structures share $\geq 95\%$ amino acid sequence identity. Afterwards,
all PDB candidates were downloaded from the Protein Data Bank, and
aligned using PyMOL using the `align PDB1, PDB2, cycles=0` command for
each PDB (PDB1) to all other identified PDBs (PDB2). For PDB structures
sharing less than 1.0 Å RMSD, the structure with the least number of
missing densities was used. In the case if any two structures shared
less than 1.0 Å RMSD and the same number of residues without missing
densities, the structure with the lowest resolution was given
preference. As a side note, for the majority of class I fusion protein
structures available within the PDB, the preponderance of determined
structures were of a single domain, most often as a monomeric domain,
which eliminated most viral fusion protein candidates from the final
dataset used to train and analyze the AxIEM model(s).\
Superimposition sometimes required reordering chain identification to
preserve residue ordering. Table 1 represents the PDB files originally
downloaded from the [Protein Data Bank](www.rcsb.org), the original
chain order, and the new chain order. If chain reordering was necessary,
the following script was used to reorder the chain IDs.\
The notation '`\backslash`' indicates that the following line should be
entered on the same line when entered in the terminal.

      python reorder_pdb_chains.py <input.pdb> <output.pdb> \
      --new_chain_order=NEW_CHAIN_ORDER \
      --new_chain_id=NEW_CHAIN_ID --preserve

After reordering the chains, the sequence of each conformation's monomer
(including either the protomer FASTA file or the two cleaved attachment
and fusion domain chains concatenated together as one FASTA file) was
aligned using Clustal Omega. All residues which aligned were considered
for epitope prediction. Any residues that were not present in all PDB
structures were removed from the respective PDB structure.

## Identification and threading of consensus sequence onto native models {#sec:msa}

Given that some of the conformations were engineered to be
conformationally stable, the original sequences of the experimentally
determined PDB ensemble did not share 100% sequence identity. Therefore,
the consensus sequence of full-length protein isolated from human hosts
was determined and threaded onto each native model backbone, that is,
the consensus sequence was used to replace the original PDB sequence so
that the amino acid identity of models were identical. In all cases the
consensus sequence shared $\geq 98\%$ sequence identity.

To generate the consensus sequence, a multiple sequence alignment was
performed using a locally-installed version of Clustal Omega
(<http://www.clustal.org>) to align all full-length sequences obtained
from the NCBI Virus, or other specialized NCBI-sponsered database ().
Sequences were were initially downloaded as a Protein FASTA file and
then aligned using the following command.\

      clustalo -i <sequence.fa> -o <sequence.aln> -t Protein --infmt=fa

Next, the consensus sequence of each multiple alignment was obtained
using `EMBOSS v.6.6.0.0` with the `cons` package
(<ftp://emboss.open-bio.org/pub/EMBOSS/>).\

      cons -sequence <sequence.aln> -outseq <sequence.cons>

::: {#tab:access_dates}
  ------------- ------------------ ------- ------- -- --
    Database      Date Accessed                       
     (taxid)                                          
     (taxid)     Collection Dates                     
    Sequences                                         
   NCBI Virus       07/13/2020                        
   \(694009\)                                         
    \(9605\)                          9               
   NCBI Virus       07/13/2020                        
   \(2697049\)                                        
    \(9605\)           all          9,339             
   NCBI Virus       09/10/2020                        
   \(208893\)                                         
    \(9605\)           all           767              
   NCBI Virus       09/10/2020                        
    \(11676\)                                         
    \(9605\)           all          80760             
   NCBI Virus       09/10/2020                        
   \(186538\)                                         
    \(9605\)           all          1700              
                                                      
    Resource        09/10/2020                        
      H3N\*           Human          all    28519     
                                                      
    Resource        09/10/2020                        
      H7N\*           Human          all     108      
  ------------- ------------------ ------- ------- -- --

  : Accession Dates of Virus Sequences. Default parameters were used
  unless noted above to query for sequences of each viral fusion
  protein.
:::

To thread the consensus sequence over each PDB, the consensus sequence
was first aligned to each of the native PDB sequences. This required
that the consensus sequence was concatenated in triplicate since the
consensus sequence only represented the full-length sequence of a viral
fusion protein monomer, while the PDB sequence represented a trimeric
sequence. Each PDB FASTA files was obtained using the following command:

      python get_fasta_from_pdb.py <pdb>

The consensus and full-length native sequences were aligned using
Clustal Omega (`<virus>_cons.fasta` and `<virus>_cons.aln`). Afterwards,
a [`grishin`
file](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/Grishan-format-alignment)
was created using the consensus sequence as the target sequence for all
native templates, and the Rosetta Partial Thread application was used to
assign coordinates to the consensus sequence. For an alternative
protocol capture on partial threading, see Section 2 of [this
tutorial](http://staging.meilerlab.org/wp-content/uploads/2022/02/rosetta_cm_tutorial.pdf).
The output PDB model was used to obtain all features, except for the
Rosetta REU residue score that requires minimization for scoring (as
described in the following section), and contains the renumbered
sequences as listed in .\

      /path/to/rosetta/main/source/bin/partial_thread.linuxgccrelease \
      -database /path/to/rosetta/main/database -in:file:fasta <sequence.cons> \
      -in:file:alignment cons_<pdb>.grishin -in:file:template_pdb <pdb>

::: {#tab:designed-positions}
                                Viral protein                                Chain   Residues considered for design
  ------------------------------------------------------------------------- ------- --------------------------------
                                EBOV Zaire GP                                  A    
                                                                               B    
                                                                               C    
                             influenza A H3 HA2                                A                387-499
                                                                               B                958-1070
                                                                               C               1529-1641
                              influenza A H7 HA                                A    
                                                                               B    
                                                                               C    
                                  HIV-1 Env                                    A    
                                                                               B    
                                                                               C    
                                    RSV F                                      A    
                                                                               B    
                                                                               C    
                                SARS-CoV-1 S                                   A    
                                                                                    
                                                                               B    
                                                                                    
                                                                               C    
                 3001-3012, 3019-3172, 3183-3318, 3345-3614                         
                                SARS-CoV-2 S                                   A    
                                                                                    
                                                                                    
                                                                               B    
                                                                                    
                                                                                    
                                                                               C    
   2746-2755, 2763-2788, 2810-2874, 2881-2987, 2988, 2995-3000, 3037-3044,          
      3049-3061, 3068-3166, 3187-3218, 3236-3356, 3360-3373, 3402-3692              

  : Residue positions considered for design. All PDB models within an
  ensemble are numbered identically, and all chain identifiers from the
  initial model are eliminated. Chain identification denote individual
  monomers of Class I fusion proteins. Residue numbering is based off of
  the threaded model (*i.e.* the $\texttt{<pdb>\_threaded.pdb}$ model),
  or rather a residue's position in the full-length consensus sequence.
  For Class I fusion proteins, only residues that are present in all
  three protomers were considered for design. Residues not present
  within all protomers were kept in the native model were allowed to
  repack (re-position) their side chains during design. HIV-1 Env was
  not subjected to threading due to low consensus sequence identity, and
  the numbering of the models used retained their original PDB
  numbering.
:::

## Energy minimization and scoring to obtain the Rosetta REU residue score {#sec:relax}

*Note*, this section requires the use of the Rosetta protein structure
prediction and design modeling suite, which is available by license for
free for non-commercial purposes, although a commercial license is
available. If you are new to using Rosetta or would like to learn more
about how install and use Rosetta, please start
[here](https://www.rosettacommons.org/docs/latest/getting_started/Getting-Started).
For the following code blocks, you will need to replace `/path/to/` with
the correct directory path to where you have installed Rosetta on your
own machine. Depending on the build you installed, you will need to
replace `linuxgccrelease` with the release version you installed.\
The threaded template models were subjected to constrained Rosetta
FastRelax to generate 50 relaxed models.

      /path/to/rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
      @relax.flags -s <threaded pdb> -scorefile <pdb>_relaxed.fasc
      
      ---------------relax.flags---------------
       -database /path/to/rosetta/main/database/
       -linmem_ig 10
       -in:file:fullatom
       -in:detect_disulf false
       -relax:fast
       -relax:constrain_relax_to_start_coords
       -out:file:fullatom
       -out:suffix _relax
       -use_input_sc
       -nstruct 50
      -----------------------------------------

The relaxed model with the combined lowest total energy score and lowest
$C_{\alpha}$ root mean square deviation (RMSD) to the threaded PDB
structure were selected as the input model to calculate per-residue REU
as follows:

      python /path/to/rosetta/tools/protein_tools/scripts/score_vs_rmsd.py  \
      -n <threaded_pdb> -c ca -t total -o <pdb>_sc_rmsd.tab <pdb>_threaded_relax_*pdb
      
      cat <pdb>_sc_rmsd.tab| tail -50 | sort -k2 -k3 | head -1 > low_model.txt
      cat low_model.txt | awk '{system("cp "$1" <pdb>_relaxed.pdb")}'
      rm low_model.txt

The per-residue energy scores were obtained using the per-residue total
energies of the Rosetta `score_jd2` output score file,
(`<pdb>_relaxed.sc`):

    /path/to/rosetta/main/source/bin/score_jd2.linuxgccrelease \
    -s <pdb>_relaxed.pdb -ignore_unrecognized_res \
    -out:file:scorefile <pdb>_relaxed.sc

# Assignment of conformation dependent epitope residues

An epitope reside is first defined here as any residue that has been
annotated as an epitope by the [Immune Epitope DataBase
(IEDB)](https://www.iedb.org), [Influenza Research
Database](https://www.fludb.org)'s Immune epitope search, or the [HIV
Molecular Immunology
Database](https://www.hiv.lanl.gov/components/sequence/HIV/featuredb/search/env_ab_search_pub.comp)
that is associated with a PDB structure. (). IEDB searches used the
filters 'Positive Assays only', 'Epitope Structure: Discontinuous', 'No
T cell assays', 'No MHC ligand assays', and 'Host: Homo sapiens
(human)'. Influenza epitope searches used the filters 'Virus Type A',
'Subtype H3 or H7', 'Protein HA, Segment 4', 'Experimentally Determined
Epitopes', 'Assay Type Category and Result B-cell Positive', and 'Host
Human'. HIV epitopes include epitopes as listed in the interactive
epitope maps as of 1 June, 2020.\
To determine each epitope residue's conformation specificity, a residue
must have at least one PDB structure of an antibody-antigen complex
where it has been annotated as an epitope residue (*i.e.* it has an IEDB
ID or is listed as an epitope on the HIV DB gp160's epitope interaction
map), and that the PDB of the annotated antibody-antigen complex, when
aligned to each monomer or chain of the benchmark protein models,
results in zero overlap of the antibody with each benchmark protein
model. Checking for overlap was performed as follows. *i*) For each PDB
antibody-antigen complex associated with an IEDB ID or HIV epitope, the
antigen and the antibody were created as independent PyMOL objects,
let's say labeled as objects `antigen` and `antibody`. *ii*) Three PyMOL
objects were created for each AxIEM benchmark PDB of a viral fusion
protein, with each object containing the residues identified to be
present in the antigen of the antibody-antigen complex, labeled as
`objA`, `objB`, and `objC`, respectively. *iii*) The `antigen` object
was first aligned to each `objX` object. *iv*) Next, the `antibody`
object was aligned to `antigen` with respect to it's aligned position to
`objX`. *v*) If no atoms of the `antibody` object came within 3 Å of any
atoms present in the AxIEM PDB model, the residues within `objX` were
considered to be a viable conformation-dependent epitope. There were
often multiple (subunit) antibody-antigen complexes associated with each
IEDB ID, and if any one of the representative complexes met the criteria
in step *v*, those residues for the given monomer/protomer were assigned
as epitope residues due to the potential of differing antibody binding
angle.\

::: {#tab:epitopes}
  ------------------------------------------- ----------------- ------------------------------------------ -----------------
                 Viral Protein                   Epitope ID                  Epitope Residues              
                 conformations                                                                             
                 Zaire EBOV GP                     442029                 N550, D552, G553, C556           
                                                   534853                       A526, I527                 
                                                   534854                                                  
                                                                                                           
                                                   534855                                                  
                                                                                                           
                                                   539006                       N550, D552                 
                                                   606556                          G528                    
                                                   857622                 N550, D552, G553, C556           
                                                   933255                    A148, G149, I532              
                                                   933256                       G118, T144                 
                                                   933257                       G149, I532                 
                                                   933258                    A525, I527, I532              
                                                   933259                    I185, I527, I532              
                                                   933260                 K115, D117, G118, T144           
                                                   933263                    R64, I527, I532               
                                                   933264                  S46, D49, G118, T144            
                                                   985426                                                  
                                                   985702                       P116, D117                 
                                                   1063108       A525, A526, I527, G528, L529, A530, W531  
               influenza H3 HA2                    189321                                                  
                                                    1HTM                                                   
                                                   580002                    Q388, I391, I394                    1HTM
                                                   580003                                                  
                                                    1HTM                                                   
                                                   742477                          I391                          1HTM
                influenza H7 HA                     H7.5                                                   
                                                    3M5G                                                   
                                                 580003^\*^                                                
                                                    6MLM                                                   
                                                   886618                                                  
                                                    6MLM                                                   
                influenza H7 HA                    886619                                                  
                                                                                                           
                                                    6MLM                                                   
                                                   886620                                                  
                                                                                                           
                                                    6MLM                                                   
                                                   952484                          G151                          6MLM
                   HIV-1 Env                       164069              C119, V120, L122, M434, P437           6U0L, 6U0N
                                                   489886               E87, N88, T90, P238, P240          
                                                   164067         C119, V120, L122, T198, 199, A200, 201      6U0L, 6U0N
                                                    16470                                                  
                                                 6U0L, 6U0N                                                
                                                   164071                                                  
                                                                                                           
                                                 6U0L, 6U0N                                                
                                                   164073                                                  
                                                                                                           
                                                                                                           
                                                   164094                                                  
                                                 6U0L, 6U0N                                                
                                                   164099                                                  
                                                 6U0L, 6U0N                                                
                                                   227937                                                  
                                                 6U0L, 6U0N                                                
                                                   534824               Q82, E83, I84, V245, Q246             6U0L, 6U0N
                                                   489875        D325, I326, R327, Q328, H330, T413, P415  
                     RSV F                         186804                                                  
                                                                                                           
                                                    3RKI                                                   
                                                    77299                          I266                    
                                                   429158              T50, L305, G307, I309, D310         
                                                   566539                    K271, L467, K470                    3RKI
                                                   566540                                                  
                                                    4MMS                                                   
                                                   581507                       N175, D263                 
                                                   581508                                                  
                                                                                                           
                                                    3RKI                                                   
                                                   581509                      S173, T174,                 
                                                   581510                    S173, T174, N175              
                                                   581511                       T174, D263                 
                                                   581512                 T174, N175, D194, D263           
                                                   591404                       G307, D310                 
                                                   606552                                                  
                                                                                                           
                                                    3RKI                                                   
                                                   912903                                                  
                                                                                                           
                                                                                                           
                  E294, E295                        3RKI                                                   
                     RSV F                         969092                                                  
                                                    3RKI                                                   
                                                   753466                                                  
                                                    4MMS                                                   
                  SARS-CoV S                        76972                          D462                     6NB7, 6NB6(B,C)
                                                    77442                                                  
                                                                                                           
                                                                                                           
                                               6NB7, 6NB6(B,C)                                             
                                                    77444                                                  
                                               6NB7, 6NB6(B,C)                                             
                                                   420672              K344, F360, Y442, L472, D480         6NB7, 6NB6(B,C)
                                                   420673                          N479                     6NB7, 6NB6(B,C)
                                                   910052                 G446, P462, D463, Y475            6NB7, 6NB6(B,C)
                                                   1074318                         D480                     6NB7, 6NB6(B,C)
                                                   1074319                K439, G446, S461, D463            6NB7, 6NB6(B,C)
                 SARS-CoV-2 S                      997006                                                  
                                                                                                           
                                                                                                           
                                                                                                           
                                                   1074327                                                 
                                                                                                           
                                                                                                           
                                                                                                           
                                                   1075135             R346, Y449, N450, L452, S494        
                                                                                                           
                                                   1075136                                                 
                                                                                                           
                                                                                                           
                                                                                                           
                                                   1083498                                                 
                                                                                                           
                                                                                                           
                                                                                                           
                                                                                                           
                                                   1087140       Y449, L492, Q493, S494, G496, Q498, Y505  
                                                                                                           
                                                                                                           
                                                   1097186                                                 
                                                                                                           
                                               all but 7CAI(B)                                             
                                                   1087266                                                 
                                                                                                           
                                                                                                           
                                                                                                           
                                                                                                           
                                                                                                           
                                                   1087267                                                 
                                                                                                           
                                                                                                           
                                                                                                           
                                                   1087269                                                 
   C379, Y380, G381, V382, S383, P384, T385,                                                               
         K386, L390, F392, D428, T430                                                                      
                 7CAI(B), 7CAK                                                                             
                                                   1087820                      D428, F429                       6VXX
                                                   1087821                         N354                          6VXX
                                                   1125015                   A372, F374, C379                    6vXX
                                                   1125016       F374, S375, T376, F377, C379, F392, D427        6VXX
                  SARS-CoV-2                       1181325             Y449, Y453, L492, Q493, S494        
                                                                                                           
                                                                                                           
                                                   1307796                                                 
                                                                                                           
                                                                                                           
                                                                                                           
                                                   1309150                                                 
                                                                                                           
                                                                                                           
                                                                                                           
                                                   1310037                      Y449, Q493                 
           7BYR(C), 7C2L(C), 7CAI(C)                                                                       
                                                   1310038                                                 
   Y453, Q493, S494, Y495, G496, T500, Y505                                                                
                6X29, 6X2B(A),                                                                             
             7BYR(B,C), 7C2L(A,B),                                                                         
                    7CAI(C)                                                                                
                                                                                                           
  ------------------------------------------- ----------------- ------------------------------------------ -----------------

  : Residues classified as epitopes. All residues listed were annotated
  to be experimentally determined epitope contacts that are present in
  all PDB models used for the AxIEM benchmark. Residue numbering refers
  only to Chain A or the original consensus sequence position number.
  The column 'Excluded Conformations' refers to any conformations (PDB
  ID) for which those residues did not meet the criteria to be
  classified as an epitope. The `AxIEM.data` file contains each
  residue's epitope label, or classifier, with a `1` indicating that
  residue was assigned to be a conformation-specific epitope residue or
  `0` if not. With the exception of HIV-1 Env, all proteins use
  contiguous numbering instead of restarting the same numbering scheme
  with each chain so that the labels had to be mapped to the correct
  position index within the dataset and models' PDB numbering schema,
  which is reflected in `AxIEM.data` but not in the table.
:::

# Description of AxIEM

## Benchmark with Discoscope and Ellipro 2.0

The annotated dataset `AxIEM.data` was constructed by concatenating the
output from Sections 1.2 and 1.3 to generate each residues REU,
CP~RMSD~, and NV per-residue features, along with the manually curated
classifier assignment with each fusion protein's list of PDBs and virus
type in the same order. Each residue in the dataset, represented by a
single row, was labeled with the virus name, PDB ID and residue number
according the models used for the benchmark, located in the `pdb_models`
directory. To compare against [Discoscope 2.0](tools.iedb.org/discotope)
and [Ellipro](tools.iedb.org/ellipro), each benchmark PDB model was
uploaded to the respective server and submitted using default
parameters. Afterwards the associated classifier label and each method's
prediction score were concatenated in the same order as listed in
`AxIEM.data` to form the benchmark datasets `Discotope.data` and
`Ellipro.data`.\
The first stage of the benchmark uses the python script
`run_AxIEM_benchmark.py` to perform three tasks: *i*) calculate the
Neighbor Sum feature as a weighted sum of all per-residue {REU,
CP~RMSD~, NV} with each feature scaled to the range \[0,1\], *ii*) run
leave-one-out tests where each viral fusion protein and its PDB
ensemble's features are excluded from the AxIEM dataset to use as the
test set for building a linear regression, Bayes Classifier, and Random
Forest Classifier model, and *iii*) analyze single model and cumulative
performance of the aforementioned models, Discotope 2.0, and Ellipro.
Given that the Neighbor Sum feature is dependent on the upper boundary
used to calculate the sigmoidal weights of surrounding residues, the
script performs leave-one-out tests and cumulative analysis for upper
boundary radii ranging from 8Å to 64Å in 8Å increments as well as using
the protein size (as number of residues).\

    python src/run_AxIEM_benchmark.py --data AxIEM.data --discotope Discotope.data \
    --ellipro Ellipro.data --individual_summary leaveoutsummary.txt \
    --roc_curves cumulativerocs.txt --summary cumulativestats.txt \
    --features features.txt

The `leavoutsummary.txt` output file contains the ROC AUC values of each
leave out test for each upper boundary condition used to generate the
associated linear regression, Bayes classifier, and Random Forest model
as well as the prediction ROC AUC values of Discotope and Ellipro for
each protein PDB ensemble, summarized in Fig2A. The `cumulativerocs.txt`
output file contains the false positive rate, true positive rate, and
associated threshold for the total performance, i.e. all leave-one-out
predictions, for each model type where model type here is defined either
as all per-residue predictions made bythe linear regression, Bayes
Classifier, and random forest classifier models built for each upper
boundary radius or as all Discotope or Ellipro per-residue predictions,
as summarized in Fig2B and Fig2C. The `cumulativestats.txt` output file
includes the ROC AUC, Youden's J statistic, and maximum Matthew's
correlation coefficient for each ROC curve in `cumulativerocs.txt`. The
`features.txt` output file contains each residue's feature set {REU,
CP~RMSD~, NV, NS~u~}, with NS~u~ being the upper boundary radius
dependent Neighbor Sum score, for plotting distributions.
