
# bacterialTSM 🧫
Bacterial metabolite profiles may contain features that are specific to particular taxonomic ranks. In Wei Chen et al., _Universal, untargeted detection of bacteria in human tissues using spatial metabolomics_, we demonstrate that ~350 markers can be identified from a set of ~600 individual isolates, and validated across various independent datasets. We propose that these taxon specific markers (TSMs) may be a suitable means of identification of the presence of certain taxons within clinical samples.

<img src="img/Species-mz-734_546-Clostridium difficile.png" alt="Title page of the 1758 edition of Linnaeus's Systema Naturæ" align="center" caption="Systema Naturæ, from Wikipedia">

>The figure above shows, on the left hand side, the intensities of m/z 734.546 divided into the 7 phylogenetic levels, and coloured according to the group in each which is most different to all others. The ROC curves for each rank are shown in the right-hand axes, with species obtaining a TPR of 0.8 and a FPR of 0. In conjunction, the plots demonstrate that this variable is highly specific to C. difficile.

## Publication
Wei Chen, et al., Universal, untargeted detection of bacteria in human tissues using spatial metabolomics, Preprint

## Data
The raw data referenced in the publication above is available in MetaboLights [2], [MTLBS10328](https://www.ebi.ac.uk/metabolights/MTBLS415). This repository contains a subset of the pre-processed data, sufficient to run the example code provided.

The example data contains a `[597 x 17]` matrix of 17 variable intensities, as well as a `[17 x 1]` vector of m/z values. The phylogeny of the 597 observations is provided as a `{597 x 7}` table, with the 7 taxonomic ranks being Gram, Phylum, Class, Order, Family, Genus and Species. Other taxonomic ranks can be included, but it is necessary that taxonomy is provided from the least to the most specific, i.e. least specific in the first column. The `VariableNames` property of the table contains the taxonomic ranks.

## Code
The code and example data in this repository explains the process of determining taxon specific markers (TSMs), i.e. assessing the specificity of a variable at each taxonomic rank. The `@TSM` class contains various methods which are each described below. The output from the class is a determination as to whether a variable is specific and at what rank.

### Usage
A subset of data can be found in `data/Example_TSM.mat`, and a workflow is presented in `Workflow.m`

### Class `@TSM`
`TSM(variableName,variableIntensities,phylogeny,alpha);`
This class contains all of the methods required to determine if a single variable is specific for a particular member of any taxonomic rank. The required inputs are the `variableIntensities` arranged as a `[n x 1]` vector, a table corresponding to the observations' `phylogeny`, along with the `variableName` and the significance level, `alpha`. The phylogenetic table should be sized `{n x p}` with `p` taxonomic ranks, set as the table's variable names. All methods are calculated automatically using the `TSM.runall` method (individual methods explained further below). Once completed, the `TSM.checkSpecificity` method applies thresholds to determine if a variable is suitably specific; the relevant results are saved in the `TSM.difftab` property.

```MATLAB
# Load the required data
load('data/Example_TSM.mat');

# Create object
t = TSM(name,x,y,alpha);

# Run all methods
t = t.runall;

# Test for specificity
t = t.checkSpecificity(tpr,fpr,numDiff,minObs)

# Output result
t.summary;

```
Class methods are as follow:

#### TSM.univariate
`TSM.univariate('anova');` For each of the taxonomic ranks, calculate the p value using ANOVA. A p value less than the alpha value is considered a significant difference, and is further analysed in the post hoc method. Results are stored in the `uv` property.

#### TSM.posthoc
`TSM.posthoc;` The p value identifies if there is a significant difference across a taxonomic rank, but not which groups differ. Using Tukey's honestly significance difference (HSD) test, the differences between each group in each taxonomic rank are assessed. Where the confidence intervals for each difference comparison do not contain 0, the individual difference for those groups is considered significant. Results are stored in the `ph` property, primarily consisting of the arrays of significant comparisons.

#### Count comparisons
`TSM.countcomparisons;` The post hoc test determines which groups are significantly different. This method first collates the count of differences for each group. The most different group within each taxonomic rank is determined (ties are settled by the group with the highest mean). At species level, if the most different species is either 'Genus sp' or 'Genus spp', these are skipped. The number of siblings for each taxon's winning group is also determined. From each taxon, the group with the highest number of significant differences is summarised in the results table, `difftab`. This table also shows the UV p values, the number of observations for the relevant group, its number of siblings and the total number of differences. 

#### Check consistency
`TSM.checkConsistency;` The table contains the name of the groups which have the most differences across each taxonomic rank. Because the univariate and post hoc methods are performed independently for each rank, the most different groups may not be taxonomically compatible, e.g. the most different species may not be a member of the most different genus. This method works from bottom to top (species up) and determines if the _species_ is a member of _genus_. Where the entry above is incomplete (because p > alpha), consistency is assumed. The consistency at the top level (gram) is not determined. Also added to the table is the MaxGroups column, which is the total number of unique entries in each taxonomic rank.

#### Determine ROC
`TSM.determineROC;` The true and false positive rates are used as markers for sensitivity and specificity of a variable. Using the most different group for each rank, the ROC curve is determined for that group against all other observations. These values are added to `difftab`.

#### Check specificity
`TSM.checkSpecificity(tpr,fpr,numDiff,minObs);`
The calculations thus far determine the most likely candidate at each taxonomic rank, but arbitrary thresholds are used to filter the list to include only the most specific/sensitive. Each of the parameters is described below.

| Parameter | Example | Description |
| :-        | :-      | :-          |
| `tpr` | `0.80` | The minimum true positive rate, determined from the ROC curve. |
| `fpr` | `0.10` | The maximum false positive rate, determined from the ROC curve. |
| `numDiff` | `[0 0 0 0 2 2 2]` | Number of differences, expressed as the MaxGroups-n for each of the `p` taxonomic ranks. A variable can be considered specific if it different to all of the other n-1 groups at that taxon. To allow some leeway at lower taxons, the corresponding value of numDiff can be increased from 0. In this example, variables should be different to `[MaxGroups-1 - 0]` for gram, phylum, class, order variables but `[MaxGroups-1 - 2]` for family, genus and species levels. Thus for species level, this is Differences >=230 (max 233) and for class level this is Differences >= 10 (max 11). The exceptions are increased for lower taxonomic ranks for situations where the lower number of observations makes it harder to differentiate between higher intensity groups that are significantly different. |
| `minObs` | `3` | The minimum number of observations; a group must have this many observations to be considered a taxon specific marker. |

#### Summary
The `TSM` class and its methods help to determine if a single variable can be considered to specific to a particular taxon. A summary for each variable can be printed using the `TSM.summary` method, and for summarising a series of variables the outputs of the `TSM.makeTable` method can be concatenated, such as for the 16 example variables provided in this repository

| Variable | Taxon| Name| NumObs| NumSiblings | Differences | AUC | FPR | TPR | p |
| :-       | :-   | :-  | :-    | :-          | :-          | :-  | :-  | :-  | :-|
| 227.2    | Phylum | Fusobacteria | 15 | 2 | 4 | 0.97948 | 0.0086957 | 0.86667 | 2.7884e-76 |
| 228.21 | Phylum | Fusobacteria | 15 | 2 | 4 | 0.97901 | 0.0069565 | 0.8 | 7.4275e-80 |
| 381.28 | Phylum | Bacteroidetes | 50 | 2 | 4 | 0.98422 | 0 | 0.82 | 4.9158e-149 |
| 609.49 | Phylum | Bacteroidetes | 50 | 2 | 4 | 0.938 | 0.0055556 | 0.86 | 9.0515e-89 |
| 610.49 | Phylum | Bacteroidetes | 50 | 2 | 4 | 0.90944 | 0.0037037 | 0.82 | 1.3639e-87 |
| 849.73 | Phylum | Bacteroidetes | 50 | 2 | 4 | 0.95637 | 0.0092593 | 0.84 | 5.6995e-80 | 
| 840.77 | Class | Bacteroidetes | 40 | 1 | 10 | 0.98743 | 0 | 0.975 | 4.7287e-57 |
| 852.74 | Class | Negativicutes | 7 | 2 | 10 | 0.92539 | 0 | 0.85714 | 4.4933e-96 |
| 608.40 | Order | Vibrionales | 3 | 6 | 22 | 1 | 0 | 1 | 1.7298e-185 |
| 648.56 | Order | Rhodospirillales | 7 | 2 | 22 | 0.92122 | 0.0068611 | 0.85714    | 3.0087e-42 |
| 572.48 | Family | Bacteroidaceae |25 | 3 | 45 | 0.98556 | 0.023009 | 0.88 | 2.3762e-90 |
| 649.50 | Family | Bacteroidaceae | 25 | 3 | 46 | 0.9635 | 0.017699 | 0.88 | 1.27e-115 |
| 560.42 | Genus | Veillonella | 4  | 1 | 79 | 0.99744 | 0.003413 | 1 | 5.2356e-19 |
| 381.34 | Species | Corynebacterium striatum | 3 | 5 | 230 | 0.9983 | 0.0017036 | 1 | 1.518e-46 |
| 603.26 | Species | Delftia acidovorans | 4 | 1 | 230 | 0.99616 | 0.0051195| 1 | 2.6342e-17 |
| 733.54 | Species | Clostridium difficile | 5 | 9  | 232 | 1 | 0 | 1  | 1.0628e-93 |
| 734.55 | Species | Clostridium difficile | 5 | 9 | 232 | 1 | 0 | 1 | 1.4195e-118 |


### Dependencies
The code was written using MATLAB 2023a.
Required toolboxes are:
- Statistics and machine learning

## References
1. Wei Chen et al. Universal, untargeted detection of bacteria in human tissues using spatial metabolomics. Preprint. 2024
2. Yurekten O, Payne T, Tejera N, et al. MetaboLights: open data repository for metabolomics. Nucleic Acids Research. 2024 Jan;52(D1):D640-D646. DOI: 10.1093/nar/gkad1045. PMID: 37971328; PMCID: PMC10767962.
