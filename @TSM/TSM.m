classdef TSM
    %TSM Determine the specificity of a variable
    %   Provide spectral intensities and the observations' phylogeny.
    %   Univariate analysis is peformed to determined differences between
    %   groups, post-hoc analysis assesses which groups are different, and
    %   subsequent stages determine a ROC curve as well as assessing the
    %   consistency of specificity across the phylogenetic tree.

    properties
        v       % name of the variable
        x       % [n x 1] array of (normalised) spectral intensities
        y       % {n x p} table of observations' phylogeny
        uv      % univariate results
        ph      % posthoc results
        difftab % table of differences.
        specific% indicator at which rank the v


    end

    properties (Hidden = true)

        tax     % taxonomic rank names
        numT    % number of ranks

        % Thresholds for specificity; see TSM.checkSpecificity
        alpha   % scalar, significance level
        fpr     % scalar between 0-1, false positive rate threshold
        tpr     % scalar between 0-1, true positive rate threshold
        numDiff % [1 x p] number of difference exceptions
        minObs  % scalar, minimum observations per taxonomic rank

    end

    methods
        function obj = TSM(var,x,y,alpha)
            %TSM Initialise object
            %   Check that inputs are of the correct size

            % Error checking
            assert(isnumeric(x),'Spectral intensities should be numeric');
            assert(isa(y,'table'),'Provide phylogeny as a table');

            % Reshape x to standard [n x 1] size
            x = reshape(x,[],1);
            sz = numel(x);

            % Expect table and data sizes to match
            assert(size(y,1) == sz,'Arrays are of incorrect size');

            % Expect alpha to be a scalar less than 1
            alfun = @(x) all(isscalar(x) & x > 0 & x < 1 & numel(x) == 1);
            assert(alfun(alpha),'Value of alpha is invalid');

            % Save properties
            obj.v = var;
            obj.x = x;
            obj.y = y;

            % Taxonomic rank names
            obj.tax = y.Properties.VariableNames;
            obj.numT = numel(obj.tax);

            % Significance level
            obj.alpha = alpha;

            % Run all of the methods here
            obj = obj.runall(obj.alpha);


        end

        function [obj] = runall(obj,alpha)
            %runall Complete all stages of the calculation

            % Use default value stored in object
            if nargin == 1
                alpha = obj.alpha;
            end

            obj = obj.univariate('anova');
            obj = obj.posthoc(alpha);
            obj = obj.countcomparisons;
            obj = obj.checkConsistency;
            obj = obj.determineROC;

        end

        function [obj] = univariate(obj,method)
            %univariate Run univariate analysis
            %   Calculate for each taxonomic group

            % Array for p values and stats structures
            pval = NaN(obj.numT,1);
            sts = cell(obj.numT,1);

            % Define the function
            switch method
                case 'anova'
                    fcn = @(x,y) anova1(x,y,'off');
                case {'kruskalwallis','kw'}
                    fcn = @(x,y) kruskalwallis(x,y,'off');
                otherwise
                    error('No other methods specified');
            end

            % Loop through each taxonomic level
            for n = 1:obj.numT
                [pval(n),~,sts{n}] = fcn(obj.x,obj.y{:,n});
            end

            % Save the results
            obj.uv.method = method;
            obj.uv.fcn = fcn;
            obj.uv.pval = pval;
            obj.uv.sts = sts;

        end

        function [obj] = posthoc(obj,alpha)
            %posthoc Run posthoc test using the results from UV analysis
            %   The output is a table showing significant differences
            %   between groups, with their associated p value.

            % This is the post-hoc analysis which uses honestly
            % significant difference (Tukey's HSD).  The outputs of
            % multcompare are: c (each row denotes a comparision betweeen
            % groups A and B, which are the first two columns. The next
            % three columns show the lower confidence limit, the
            % difference in means, and the higher confidence limit.
            % If the confidence limit contains zero, the difference of
            % this comparison is not significant); m (group mean values
            % and standard errors); gn (group names).

            % Alpha
            if nargin == 2
                obj.alpha = alpha;
            end
            alpha = obj.alpha;

            % Temporary array for storing each taxon's results
            tmp = cell(obj.numT,1);

            % For each taxonomic level
            for n = 1:obj.numT

                % Skip insignificant differences
                if obj.uv.pval(n) > alpha
                    continue;
                end

                % Run posthoc
                [c,~,~,gn] = multcompare(obj.uv.sts{n},...
                    'ctype','hsd',...
                    'alpha',alpha,...
                    'display','off');

                % Expect names to match
                assert(all(strcmp(gn,obj.uv.sts{n}.gnames)),...
                    'Group names from univariate -> posthoc do not match');

                % Determine any significant differences
                ll = c(:,3) <= 0;
                hh = c(:,5) >= 0;
                ci = sum([ll hh],2) == 1;

                % Make a the list of group differences by their names
                c2 = c(ci,:);
                tmp{n} = c2(:,[1 2 6]);

            end

            % Save the results
            obj.ph.method = 'hsd';
            obj.ph.alpha = alpha;
            obj.ph.df = tmp;

        end

        function [obj] = countcomparisons(obj)
            %countcomparisons Summarise the results from posthoc analysis
            %   Determine groups which differ to most other groups

            % Summarise the findings in a single table
            c = table(cell(obj.numT,1),NaN(obj.numT,1),...
                NaN(obj.numT,1),NaN(obj.numT,1),NaN(obj.numT,1),...
                'VariableNames',{'Name','p','NumObs','NumSiblings',...
                'Differences'},...
                'RowNames',obj.tax);

            % For each taxonomic rank
            for n = 1:obj.numT

                % Is the table empty?
                if isempty(obj.ph.df{n})
                    continue
                end

                % Here run the function
                [cnt] = obj.comparisonCountUp(obj.ph.df{n},obj.uv.sts{n});

                % If the different group is `Genus sp` then we don't
                % consider as a separate species, so can skip this one here
                if n == obj.numT && ...
                        (endsWith(cellstr(cnt{1,'Label'}),' sp') || ...
                        endsWith(cellstr(cnt{1,'Label'}),' spp'))
                    continue;
                end

                % Save most different to the summary table
                c{n,'Name'} = cellstr(cnt{1,'Label'});

                % Number of observations in this group
                fx = strcmp(obj.y{:,n},c{n,'Name'});
                c{n,'NumObs'} = sum(fx);

                % How many siblings does it have?
                if n > 1
                    parent = unique(obj.y{fx,n-1});
                    fy = strcmp(obj.y{:,n-1},parent);
                    siblings = unique(obj.y{fy,n});

                    % What about `Genus sp` species? These are not
                    % considered to be species level, so should be removed
                    % from the species level sum
                    if n == obj.numT

                        % Am I a `Genus sp` species?
                        iamSp = endsWith(c{n,'Name'},' sp') | ...
                            endsWith(c{n,'Name'},' spp');

                        % Do I have any `Genus sp` siblings? (Me included)
                        childrenSp = endsWith(siblings,' sp') | ...
                            endsWith(siblings,' spp');
                        numSibSp = sum(childrenSp) - double(iamSp);
                    else
                        numSibSp = 0;
                    end

                    % Final determination
                    numSib = numel(siblings) - 1 - numSibSp;

                else
                    % Set to 1 at the top level
                    numSib = 1;
                end
                c{n,'NumSiblings'} = numSib;


                % Number of significant differences
                c(n,'Differences') = cnt(1,'Count');
            end

            % Save this as the table of differences
            c.p = obj.uv.pval;
            obj.difftab = c;

        end

        function [obj] = checkConsistency(obj)
            %checkConsistency Is the variable consistent up the tree?
            %   Starting from the bottom of the tree, determine if it
            %   is a member of the one above it, e.g. is [Clostridiales]
            %   a member of [Clostridia]?

            % Determine the tree, by reducing the size of the taxonomic
            % labels
            [~,idx,~] = unique(obj.y(:,end));
            tree = obj.y(idx,:);
            size(tree);

            consistency = NaN(obj.numT,1);

            % Maximum number of groups in each taxonomic level
            maxGroups = NaN(obj.numT,1);
            for g = 1:obj.numT
                maxGroups(g) = numel(unique(tree(:,g).Variables));
            end
            obj.difftab.MaxGroups = maxGroups;

            % What if there are no taxon entries? (i.e. completely
            % non-specific)
            if all(cellfun(@isempty,obj.difftab{:,'Name'}))
                obj.difftab.Consistency = false(obj.numT,1);
                return
            end

            % From bottom to top... except top layer
            for g = obj.numT:-1:2

                % Name of most specific group at this taxonomic level
                taxThis = obj.difftab{ g ,'Name'};
                taxPrev = obj.difftab{g-1,'Name'};

                % If the taxon above is not listed, then consistency is not
                % a problem because there were no significant differences
                % in the group above
                if isempty(taxPrev{1})
                    if isempty(taxThis{1})
                        consistency(g) = NaN;
                    else
                        consistency(g) = 1;
                    end

                    continue;
                end

                % Skip if empty...
                if ~isempty(taxThis{1})

                    % Is [Species] a member of [Genus]?

                    % Find species in the tree
                    thisLevelInTree = strcmp(tree(:,g).Variables,taxThis);

                    % Find supposed parent in the tree
                    nextLevelInTree = strcmp(tree(:,g-1).Variables,taxPrev);

                    % Is there a match?
                    if any(thisLevelInTree & nextLevelInTree)
                        consistency(g) = 1;
                    else
                        consistency(g) = 0;
                    end
                end
            end

            % Save to the table
            obj.difftab.Consistency = consistency;

        end

        function [obj] = determineROC(obj)
            %determineROC Calculate ROC curve and associated metrics
            %   ~

            % Store AUC, FPR, TPR
            obj.difftab.AUC = NaN(obj.numT,1);
            obj.difftab.FPR = NaN(obj.numT,1);
            obj.difftab.TPR = NaN(obj.numT,1);

            % Loop through each taxon; skip non consistent levels, or those
            % where no name is provided (should overlap)
            for n = 1:obj.numT

                if obj.difftab.Consistency(n) ~= 0 ...
                        && ~isempty(obj.difftab.Name{n})

                    % Determine AUC, FPR, TPR
                    [~,~,~,auc,optauc] = perfcurve(obj.y{:,n},...
                        obj.x,...
                        obj.difftab{n,'Name'});

                    obj.difftab.AUC(n) = auc;
                    obj.difftab.FPR(n) = optauc(1);
                    obj.difftab.TPR(n) = optauc(2);
                end

            end

        end

        function [obj] = checkSpecificity(obj,tpr,fpr,numDiff,minObs)
            %checkSpecificity Determine if a variable is specific
            %   Thresholds are used to determine if a variable passes the
            %   various tests, and a variable is deemed specific if all are
            %   passed.

            % Save the thresholds to the object
            obj.tpr = tpr;
            obj.fpr = fpr;
            obj.numDiff = numDiff;
            obj.minObs = minObs;

            % Check p values from each level
            chkP = obj.uv.pval <= obj.alpha;

            % Specific enough to pass the tpr/fpr tests
            chkROC = obj.difftab.TPR >= tpr & ...
                obj.difftab.FPR <= fpr;

            % Different to at least a large number of groups in each taxon
            chkGroups = obj.difftab.Differences >= ...
                obj.difftab.MaxGroups-1-reshape(numDiff,[],1);

            % At least minObs observations for this taxon
            chkObs = obj.difftab.NumObs >= minObs;

            % If there is only one group in any one taxon, we won't
            % consider the variable to be specific at this level, such as
            % only one species in a genus, or one family in an order. In
            % such cases, the variable can be considered specific in the
            % less specific taxon
            chkSib = obj.difftab.NumSiblings > 0;

            % Check consistency; only go as far down the list as the first
            % non-consistent entry, starting from the top.
            tmp = find([obj.difftab.Consistency;0] == 0,1,'first') - 1;
            chkCon = false(obj.numT,1);
            chkCon(1:tmp) = true;

            % Finally combine all metrics
            chkAll = all([chkP chkROC chkGroups chkObs chkSib chkCon],2);

            % Add to the table
            obj.difftab.Specific = chkAll;

            % Use the one at the most specific level
            if any(chkAll)
                obj.specific = find(chkAll,1,'last');
            end

        end

        function summary(obj)
            %summary Write a summary of the variable
            %   Say if the variable is specific, at what level, and the
            %   vital statistics.

            if isnan(obj.specific)

                % Variable is not specific
                txt = sprintf('Variable %0.2f is not specific at any level.',...
                    obj.variable);

            else

                n = obj.specific;

                % The variable is specific etc
                txt = sprintf(['Variable %0.2f is specific to %s (%s '...
                    'level) with a TPR of %0.2f and FPR of %0.2f. The '...
                    'variable is different to %d (out of %d other) groups ' ...
                    'at the %s level. There are %d observations for '...
                    'this group.'],...
                    obj.v,obj.difftab{n,'Name'}{1},...
                    obj.tax{n},obj.difftab{n,'TPR'},...
                    obj.difftab{n,'FPR'},obj.difftab{n,'Differences'},...
                    obj.difftab{n,'MaxGroups'}-1,obj.tax{n},...
                    obj.difftab{n,'NumObs'});

            end

            % Print out to screen
            disp(txt);

        end

        function [o] = makeTable(obj)
            %makeTable Output a single-row table
            %   Can later combine with other variables' tables

            n = obj.specific;

            % Variable names
            varNames = {'Variable','Taxon','Name','NumObs',...
                'NumSiblings','Differences','AUC','FPR','TPR','p'};

            % Either / or
            if isnan(n)
                % Change the table for non-specific variables
                o = table(obj.v,{''},{''},...
                    NaN,NaN,NaN,...
                    NaN,NaN,NaN,NaN,...
                    'VariableNames',varNames);
            else

                d = obj.difftab(n,:);

                o = table(obj.v,obj.tax(n),d.Name,...
                    d.NumObs,d.NumSiblings,d.Differences,...
                    d.AUC,d.FPR,d.TPR,obj.uv.pval(n),...
                    'VariableNames',varNames);

            end



        end

    end

    methods (Static = true)

        function [cnt] = comparisonCountUp(tab,sts)
            %comparisonCountUp Determine frequently different groups
            %   Input is a table from the posthoc function and the output
            %   of the stats structure from ANOVA

            % Convert entries from numerical values to taxon names
            l1 = sts.gnames(tab(:,1));
            l2 = sts.gnames(tab(:,2));

            % Summarise the entries
            cnt = countlabels(cat(1,l1,l2));

            % Append the mean of each group to the table, which helps to
            % break ties.
            [~,b] = ismember(cellstr(cnt{:,'Label'}),sts.gnames);
            cnt.Mean = sts.means(b)';

            % Sort the table (most differences, then largest mean)
            cnt = sortrows(cnt,{'Count','Mean'},'descend');

        end


    end

end