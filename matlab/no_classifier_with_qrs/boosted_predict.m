function y = boosted_predict(X,clf,output_type,max_round)
% FUNCTION y = boosted_predict(X,clf,output_type)
%   *** INPUTS ***
%   X:
%   clf:
%   output_type:
%
%   *** OUTPUTS ***
%   y:

if ~exist('output_type','var'); output_type = ''; end;
if ~exist('max_round','var'); max_round = []; end;

if ~numel(max_round); max_round = inf; end;

weak_learners = clf.weak_learners;
stacking = clf.stacking;

[n,p] = size(X);
X = [X,ones(n,1)];
p = p + 1;

if strcmp(clf.missingDataHandler.method,'abstain')
    % Then we have nothing to do...
elseif strcmp(clf.missingDataHandler.method,'mean')
    % Use the imputed values
    for j=1:p
        X(isnan(X(:,j)),j) = clf.missingDataHandler.imputedVals(j);
    end
elseif strcmp(clf.missingDataHandler.method,'multivariate_normal_model')
    % Then use the multivariate normal to impute values on each trial
    % First, get all of the distinct missing patterns
    [MP,~,pattern_inds] = unique(~isnan(X(:,1:end-1)),'rows');
    [pattern_index,pattern_inds] = sort(pattern_inds,'ascend');
    j = 0;
    while j < n
        fprintf('%d\n',j);pause(1e-5);
        j = j + 1;
        jStart = j;
        index = pattern_index(j);
        while j <= n && pattern_index(j) == index
            j = j + 1;
        end
        j = j - 1;
        
        inds = pattern_inds(jStart:j);
        
        exists_locs = find(MP(index,:));
        missing_locs = setdiff(1:p-1,exists_locs);
        
        tmpCov = clf.missingDataHandler.covMat;
        tmpCov = tmpCov([missing_locs,exists_locs],:);
        tmpCov = tmpCov(:,[missing_locs,exists_locs]);
        
        S12 = tmpCov(1:numel(missing_locs),(numel(missing_locs)+1):end);
        S22 = tmpCov((numel(missing_locs)+1):end,(numel(missing_locs)+1):end);
        
        impute_vals = repmat(clf.missingDataHandler.means(missing_locs),1,numel(inds)) + S12*(S22\(X(inds,exists_locs)' - repmat(clf.missingDataHandler.means(exists_locs),1,numel(inds))));
        
        X(inds,missing_locs) = impute_vals';
    end
else
    error('Unknown missingDataHandler method');
end

XM = double(~isnan(X));
yfull = zeros(n,p);yfullind=0;%min(max_round,numel(decision_stumps)));yfullind = 0;
y = zeros(n,1);
for j=1:p
    if ~numel(weak_learners{j}); continue; end;
    
    locs = find(XM(:,j));
    ycurr = zeros(n,1);
    for k=1:numel(weak_learners{j})
        ds = weak_learners{j}{k};
        if ds.boosting_round > max_round
            break;
        end
        if isfield(ds,'bias')
            ycurr(locs) = ycurr(locs) + ds.bias;
        end
        
        switch ds.type
            case 'constant'
                % There's nothing left to do
            case 'logistic'
                ycurr(locs) = ycurr(locs) + ds.alpha*(2./(1+exp(-(ds.a*X(locs,j)+ds.b)))-1);
            case 'stump'
                % First assume it's positive and then double correct the negatives
                ycurr(locs) = ycurr(locs) + ds.alpha;
                neg_locs = find(X(locs,j) < ds.threshold);
                ycurr(locs(neg_locs)) = ycurr(locs(neg_locs)) - 2*ds.alpha;
            otherwise
                error('Unknown weak learner type');
        end
    end
    
    if isfield(stacking,'mf_coeffs')
        c = stacking.mf_coeffs{j}(1) + XM*stacking.mf_coeffs{j}(2:end)';
        ycurr = ycurr.*c;
    end
    
    yfullind = yfullind + 1;
    yfull(:,yfullind) = ycurr;
    
    y = y + ycurr;
end

if strcmp(output_type,'')
    % Then don't do anything
elseif strcmp(output_type,'probability')
    y = 1./(1+exp(-y));
else
    error('Unknown output_type');
end
end