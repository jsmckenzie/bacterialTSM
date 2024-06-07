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
Univariate analysis at each taxonomic rank, using ANOVA. A p value for each rank is determined, i.e. is one of the groups significantly different from any of the others?

#### TSM.posthoc
At each taxonomic rank, if the univariate p value is less than `0.05`, then Tukey's honestly significant difference (HSD) test is performed. For each rank, the output is a count of the number of significant differences between a group and all others. If the intensities of one member of a rank (e.g. `x` in class `y`) are significantly higher than all of the (e.g.) 10 other classes, it will differ to 9 of them. Ranks that have lower intensities than other members are not considered.

#### Count differences


#### Determine ROC


#### Summarise



### Dependencies
The code was written using MATLAB 2023a.
Required toolboxes are:
- Statistics and machine learning
- ?

## References
1. Wei Chen et al.
2. Yurekten O, Payne T, Tejera N, et al. MetaboLights: open data repository for metabolomics. Nucleic Acids Research. 2024 Jan;52(D1):D640-D646. DOI: 10.1093/nar/gkad1045. PMID: 37971328; PMCID: PMC10767962.
