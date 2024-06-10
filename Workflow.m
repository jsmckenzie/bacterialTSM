%% bacterialTSM
% Workflow for the determination of taxon specific markers.

%% Load the data

clearvars;

load('data/Example_TSM.mat');


%% Create the object
% The input arguments are:
% - variable name/description, e.g. m/z value
% - [n x 1] vector of intensities
% - {n x p} table containing taxonomic information
% - [n] scalar, the significance level alpha

% Create the object, using m/z = 733.54 as an example, which is specific 
% at the species level to C. difficile
t = TSM(varName(16),vars(:,16),phyTab,alpha);


%% Methods
% The following methods need to be completed in order to determine the
% specificity of the variable:
%   1. TSM.univariate
%   2. TSM.posthoc
%   3. TSM.countcomparisons
%   4. TSM.checkConsistency
%   5. TSM.determineROC
%   6. TSM.checkSpecificity


%% Univariate
% For each of the taxonomic ranks, calculate the p value using ANOVA. A p
% value less than the alpha value is considered a significant difference,
% and is further analysed in the post hoc method. Results are stored in the
% `uv` property.

t = t.univariate('anova');


%% Post hoc using HSD
% The p value identifies if there is a significant difference across a
% taxonomic rank, but not which groups differ. Using Tukey's honestly
% significance difference test, the differences between each group in each
% taxonomic rank are assessed. Where the confidence intervals for each
% difference comparison do not contain 0, the individual difference for
% those groups is considered significant. Results are stored in the `ph`
% property, primarily consisting of the arrays of significant comparisons.

t = t.posthoc;


%% Count comparisons
% The post hoc test determines which groups are significantly different.
% This method first collates the count of differences for each group. The
% most different group within each taxonomic rank is determined (ties are
% settled by the group with the highest mean). At species level, if the 
% most different species is either 'Genus sp' or 'Genus spp', these are
% skipped. The number of siblings for each taxon's winning group is
% also determined. From each taxon, the group with the highest number of
% significant differences is summarised in the results table, `difftab`.
% This table also shows the UV p values, the number of observations for the
% relevant group, its number of siblings and the total number of
% differences. 

t = t.countcomparisons;

% Display the table
disp(t.difftab);


%% Check consistency
% The table contains the name of the groups which have the most differences
% across each taxonomic rank. Because the univariate and post hoc methods
% are performed independently for each rank, the most different groups may
% not be taxonomically compatible, e.g. the most different species may not 
% be a member of the most different genus. This method works from bottom to
% top (species up) and determines if the [species] is a member of
% [genus]. Where the entry above is incomplete (because p > alpha), 
% consistency is assumed. The consistency at the top level (gram) is not
% determined. Also added to the table is the MaxGroups column, which is the
% total number of unique entries in each taxonomic rank.

t = t.checkConsistency;

% Display the table
disp(t.difftab);


%% Determine ROC curve
% The true and false positive rates are used as markers for sensitivity and
% specificity of a variable. Using the most different group for each rank,
% the ROC curve is determined for that group against all other
% observations. These values are added to `difftab`.

t = t.determineROC;

% Display the table
disp(t.difftab);



%% Check specificity
% The calculations thus far determine the most likely candidate at each
% taxonomic rank, but arbitrary thresholds are used to filter the list to
% include only the most specific/sensitive. Each of the thresholds is
% described below.

% The minimum true positive rate
tpr = 0.80;

% The maximum false negative rate
fpr = 0.10;

% Number of differences, expressed as the MaxGroups-n for each taxonomic
% rank. A variable can be considered specific if it different to all of
% the other n-1 groups at that taxon. To allow some leeway at lower taxons,
% the corresponding value of numDiff can be increased from 0. In this
% example, variables should be different to [MaxGroups-1 - 0] for gram,
% phylum, class, order variables but [MaxGroups-1 - 2] for family, genus
% and species levels. Thus for species level, this is Differences >=230 
% (max 233) and for class level this is Differences >= 10 (max 11). The
% exceptions are increased for lower taxonomic ranks for situations where
% the lower number of observations makes it harder to differentiate between
% higher intensity groups that are significantly different.
numDiff = [0 0 0 0 2 2 2];

% The minimum number of observations; a group must have this many
% observations to be considered a TSM.
minObs = 3;


% Run the method, which checks the following at each taxonomic rank:
% - p value < alpha
% - true positive rate greater than or equal threshold
% - false negative rate less than or equal to threshold
% - different to at least the minimum number of other groups in taxon
% - at least the minumum number of observations
% - at least one sibling, e.g. more than one species in a genus
% - taxonomic consistency from the top, working downwards.
% The lowest rank that satisfies all of these conditions marks the variable
% as being a TSM, and the index is stored in `specific` property of the
% class.
t = t.checkSpecificity(tpr,fpr,numDiff,minObs);

% Display the table
disp(t.difftab);





%% Phylum level @ 609.49
% Bacteroidetes

% Create object, and run all
i = 4;
t = TSM(varName(i),vars(:,i),phyTab,alpha);
t = t.checkSpecificity(tpr,fpr,numDiff,minObs);

% Show the table
disp(t.difftab);

% Summary
% For all levels, the results from ANOVA denote a signficant difference in
% group means, with all p < alpha. At species level, the variable cannot be
% considered for the following: insufficient observations; not different to
% at least 230 other species groups; true positive rate too low. At genus
% and order levels there are no siblings (i.e. only one genus of the
% Bacteroidaceae family). Phylum level is the only taxonomic rank that
% satisfies all of the criteria.


%% Class level @ 852.74
% Negativicutes

% Create object, and run all
i = 8;
t = TSM(varName(i),vars(:,i),phyTab,alpha);
t = t.checkSpecificity(tpr,fpr,numDiff,minObs);

% Show the table
disp(t.difftab);

% Summary
% From class down to species, the phylogeny is consistent, yet the low
% number of observations from family-species prevent this variable being
% marked as specific at lower levels, despite high ROC curve values and
% differences to the maximum number of other groups. As such, the variable
% is limited to a class-level one.


%% Order level @ 648.56
% Rhodospirillales

% Create object, and run all
i = 10;
t = TSM(varName(i),vars(:,i),phyTab,alpha);
t = t.checkSpecificity(tpr,fpr,numDiff,minObs);

% Show the table
disp(t.difftab);

% Summary
% Significant p values are calculated from class down to species. The only
% species of the Roseomonas genus is 'Roseomonas sp', which is henceforth
% excluded. Also, family and genus are prevent from being specific due to
% having no siblings from the parent taxon which leaves only class and
% order, but the former is excluded due to a low TPR. Consequently, the
% variable is specific at the order rank.


%% Family level @ 649.51
% Bacteroidaceae

% Create object, and run all
i = 12;
t = TSM(varName(i),vars(:,i),phyTab,alpha);
t = t.checkSpecificity(tpr,fpr,numDiff,minObs);

% Show the table
disp(t.difftab);

% Summary
% The prevotella genus is observed to be the most different, yet this is
% not consistent with either of the groups found at family and species
% levels; as such, genus and species are excluded. At family level, a
% maximum of 46/47 differences are observed for Bacteroidaceae, and equally
% for the Bacteroidetes class. As family is lower down the tree than class,
% the variable is marked as being family specific.


%% Genus level @ 560.42
% Veillonella

% Create object, and run all
i = 13;
t = TSM(varName(i),vars(:,i),phyTab,alpha);
t = t.checkSpecificity(tpr,fpr,numDiff,minObs);

% Show the table
disp(t.difftab);

% Summary
% The lower 5 ranks have significant p values, with all also passing the
% required number of differences. At species level, this rank is exluded
% because the number of observations is too low and it is not consistent
% to the genus. Both genus and family levels pass the required tests.


%% Species level @ 733.54
% Clostridium difficile

% Create object, and run all
i = 16;
t = TSM(varName(i),vars(:,i),phyTab,alpha);
t = t.checkSpecificity(tpr,fpr,numDiff,minObs);

% Show the table
disp(t.difftab);

% Summary
% All p values are significant, and the groups with the largest number of
% differences are all parent taxons of C. difficile, which is significantly
% different to all ofhter 232 species. The FPR and TPR are 0 and 1, and
% with 5 observations, the variable is marked specific at species level.



