<img src="https://upload.wikimedia.org/wikipedia/commons/a/ae/Linnaeus1758-title-page.jpg" alt="Title page of the 1758 edition of Linnaeus's Systema Naturæ" height="300" align="right" caption="Systema Naturæ, from Wikipedia">

# bacterialTSM 🧫
> Succinct summary here

## Publication
Wei Chen, et al., Universal, untargeted detection of bacteria in human tissues using spatial metabolomics, Preprint

## Data
The raw data referenced in the publication above is available in MetaboLights [2], [MTLBS10328](https://www.ebi.ac.uk/metabolights/MTBLS415).

## Code
The code and example data in this repository explains the process of determining taxon specific markers (TSMs), i.e. assessing the specificity of a variable at each taxonomic rank. The `@TSM` class contains various methods which are each described below. The output from the class is a determination as to whether a variable is specific and at what rank.

### Usage
A subset of data can be found in `data/Example.mat`, and a workflow is presented in `Workflow.m`

### Class `@TSM`
Class overview here.
```
# Load the required data

# Define TSM thresholds
minObs = []
fpr = 0.05;
tpr = 0.80;
numDiffs = [0 0 0 0 0 1 2];

# Create object
t = TSM(name,x,y);

# Run all methods
t = t.runall;

# Test for specificity
t = t.checkSpecific(~,~,~,~)

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

| Parameter | Default | Description |
| :-        | :-      | :-          |
| `tpr` | `0.80` | The minimum true positive rate, determined from the ROC curve. |
| `fpr` | `0.10` | The maximum false positive rate, determined from the ROC curve. |
| `numDiff` | `[0 0 0 0 2 2 2]` | Number of differences, expressed as the MaxGroups-n for each taxonomic rank. A variable can be considered specific if it different to all of the other n-1 groups at that taxon. To allow some leeway at lower taxons, the corresponding value of numDiff can be increased from 0. In this example, variables should be different to `[MaxGroups-1 - 0]` for gram, phylum, class, order variables but `[MaxGroups-1 - 2]` for family, genus and species levels. Thus for species level, this is Differences >=230 (max 233) and for class level this is Differences >= 10 (max 11). The exceptions are increased for lower taxonomic ranks for situations where the lower number of observations makes it harder to differentiate between higher intensity groups that are significantly different. |
| `minObs` | `3` | The minimum number of observations; a group must have this many observations to be considered a taxon specific marker. |


#### Summarise



### Dependencies
The code was written using MATLAB 2023a.
Required toolboxes are:
- Statistics and machine learning
- ?

## References
1. Wei Chen et al.
2. Yurekten O, Payne T, Tejera N, et al. MetaboLights: open data repository for metabolomics. Nucleic Acids Research. 2024 Jan;52(D1):D640-D646. DOI: 10.1093/nar/gkad1045. PMID: 37971328; PMCID: PMC10767962.
