function [res, diagOut] = conf_index_adjustv2(res,i,CATMOS,predVT,predNT,predEPA,predGHS,predLD50,varargin)
%CATMOS_ENDPOINT_CONF_COVERAGE_ADJUST_V2 Pre-WoE endpoint confidence coverage adjustment.
%
% Purpose
% -------
% Adjust the five single-endpoint confidence indices BEFORE woe_corr runs:
%   res.Conf_index_VT
%   res.Conf_index_NT
%   res.Conf_index_EPA
%   res.Conf_index_GHS
%   res.Conf_index_LD50
%
% The adjustment is based only on the availability of experimental values
% among each endpoint model's own nearest neighbors. This is intentionally
% mild because these endpoint confidences can affect the final WoE bin.
%
% Default behavior with BonusThreshold = 3, KTarget = 5:
%   0 experimental neighbors -> penalty
%   1 experimental neighbor  -> penalty
%   2 experimental neighbors -> small penalty
%   3 experimental neighbors -> small bonus
%   4 experimental neighbors -> moderate bonus
%   5 experimental neighbors -> full bonus
%
% Example call, immediately before woe_corr:
%   [res,endpointConfCovDiag] = catmos_endpoint_conf_coverage_adjust( ...
%       res,i,CATMOS,predVT,predNT,predEPA,predGHS,predLD50, ...
%       'BonusThreshold',3,'BonusMax',0.06,'PenaltyMax',0.08);
%
% Notes
% -----
% * This helper clips endpoint confidence values to [0,1] before adjustment.
% * It preserves the original confidence values in res.Conf_index_*_preCov.
% * It does not change AD_index_* or AD_*.
% * It should be called after the initial Conf_index_* formulas are computed
%   and before woe_corr/woe_corrV76 is called.

p = inputParser;
p.addParameter('KTarget',5,@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('BonusThreshold',3,@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('BonusMax',0.1,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('PenaltyMax',0.1,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('EffectPower',1.0,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('EffectMode','additive',@(x) ischar(x) || isstring(x));
p.addParameter('BlendStrength',0.25,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('ClipOutput',true,@(x) islogical(x) && isscalar(x));
p.addParameter('StoreDiagnostics',false,@(x) islogical(x) && isscalar(x));
p.addParameter('DoNotBoostWhenADIndexZero',true,@(x) islogical(x) && isscalar(x));
p.addParameter('DoNotBoostWhenADGlobalZero',false,@(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

modelNames = {'VT','NT','EPA','GHS','LD50'};
predModels = {predVT,predNT,predEPA,predGHS,predLD50};

nModels = numel(modelNames);
diagOut = struct();
diagOut.modelNames = modelNames;
diagOut.nExp = NaN(1,nModels);
diagOut.coverageRatio = NaN(1,nModels);
diagOut.firstNeighborHasExp = false(1,nModels);
diagOut.coverageEffect = NaN(1,nModels);
diagOut.coverageScore = NaN(1,nModels);
diagOut.rawConf = NaN(1,nModels);
diagOut.clippedConf = NaN(1,nModels);
diagOut.adjustedConf = NaN(1,nModels);
diagOut.effectMode = char(opt.EffectMode);
diagOut.blendStrength = opt.BlendStrength;
diagOut.AD_index = NaN(1,nModels);
diagOut.AD_global = NaN(1,nModels);

for m = 1:nModels
    endpoint = modelNames{m};
    Pm = predModels{m};

    confField = ['Conf_index_' endpoint];
    adIndexField = ['AD_index_' endpoint];
    adField = ['AD_' endpoint];

    rawConf = get_res_scalar(res,confField,i,NaN);
    baseConf = clip01(rawConf);

    ADi = get_res_scalar(res,adIndexField,i,NaN);
    ADg = get_res_scalar(res,adField,i,NaN);

    ids = [];
    if isfield(Pm,'neighbors') && ~isempty(Pm.neighbors)
        ids = get_row(Pm.neighbors,i);
    end

    expVec = get_endpoint_exp_vector(CATMOS,endpoint);
    [nExp,coverageRatio,firstHasExp] = count_experimental_neighbors(expVec,ids,opt.KTarget);

    effect = coverage_effect(nExp,opt.KTarget,opt.BonusThreshold,opt.BonusMax,opt.PenaltyMax,opt.EffectPower);
    coverageScore = coverage_score(nExp,opt.KTarget,opt.BonusThreshold,opt.EffectPower);

    switch lower(char(opt.EffectMode))
        case 'additive'
            % Original V1 behavior: small signed offset around the existing endpoint confidence.
            adjConf = baseConf + effect;

        case 'blend'
            % Stronger behavior: blend endpoint confidence toward an availability score.
            % This can affect WoE more than additive mode because low coverage can pull
            % high endpoint confidence down substantially, and high coverage can support it.
            adjConf = (1-opt.BlendStrength)*baseConf + opt.BlendStrength*coverageScore;

        case 'hybrid'
            % Additive adjustment plus a mild blend toward the availability score.
            additiveConf = baseConf + effect;
            adjConf = (1-opt.BlendStrength)*additiveConf + opt.BlendStrength*coverageScore;

        otherwise
            error('EffectMode must be additive, blend, or hybrid.');
    end

    % Preserve the original behavior when AD_index is zero: no confidence.
    if opt.DoNotBoostWhenADIndexZero && isfinite(ADi) && ADi <= 0
        adjConf = 0;
    end

    % Optional: do not boost an endpoint that is globally outside AD.
    if opt.DoNotBoostWhenADGlobalZero && isfinite(ADg) && ADg <= 0 && adjConf > baseConf
        adjConf = baseConf;
    end

    if opt.ClipOutput
        adjConf = clip01(adjConf);
    end

    % Write adjusted confidence back into the endpoint field used by woe_corr.
    res.(confField)(i,1) = adjConf;

    if opt.StoreDiagnostics
        res = store_diag_scalar(res,['Conf_index_' endpoint '_preCov'],i,rawConf,confField);
        res = store_diag_scalar(res,['Conf_index_' endpoint '_clipCov'],i,baseConf,confField);
        res = store_diag_scalar(res,['Conf_index_' endpoint '_covAdjusted'],i,adjConf,confField);
        res = store_diag_scalar(res,['NN_exp_count_' endpoint],i,nExp,confField);
        res = store_diag_scalar(res,['NN_exp_coverage_ratio_' endpoint],i,coverageRatio,confField);
        res = store_diag_scalar(res,['NN_exp_first_has_exp_' endpoint],i,double(firstHasExp),confField);
        res = store_diag_scalar(res,['Conf_index_' endpoint '_coverage_effect'],i,effect,confField);
        res = store_diag_scalar(res,['Conf_index_' endpoint '_coverage_score'],i,coverageScore,confField);
    end

    diagOut.nExp(m) = nExp;
    diagOut.coverageRatio(m) = coverageRatio;
    diagOut.firstNeighborHasExp(m) = firstHasExp;
    diagOut.coverageEffect(m) = effect;
    diagOut.coverageScore(m) = coverageScore;
    diagOut.rawConf(m) = rawConf;
    diagOut.clippedConf(m) = baseConf;
    diagOut.adjustedConf(m) = adjConf;
    diagOut.AD_index(m) = ADi;
    diagOut.AD_global(m) = ADg;
end

end

%% ------------------------------------------------------------------------
% Helper functions
%% ------------------------------------------------------------------------

function effect = coverage_effect(nExp,KTarget,bonusThreshold,bonusMax,penaltyMax,powerVal)
% Availability-only bonus/penalty.
% For KTarget=5 and bonusThreshold=3:
%   nExp 0,1,2 => penalty
%   nExp 3,4,5 => bonus

KTarget = max(1,round(KTarget));
bonusThreshold = min(max(round(bonusThreshold),1),KTarget);
nExp = min(max(round(nExp),0),KTarget);

if nExp >= bonusThreshold
    % Scales from small bonus at threshold to full bonus at KTarget.
    denom = max(KTarget - bonusThreshold + 1,1);
    frac = (nExp - bonusThreshold + 1) / denom;
    effect = bonusMax * frac.^powerVal;
else
    % Scales from full penalty at 0 to small penalty just below threshold.
    denom = max(bonusThreshold,1);
    frac = (bonusThreshold - nExp) / denom;
    effect = -penaltyMax * frac.^powerVal;
end

end

function score = coverage_score(nExp,KTarget,bonusThreshold,powerVal)
%COVERAGE_SCORE Availability score in [0,1] for blend/hybrid modes.
% For KTarget=5, bonusThreshold=3, power=1:
%   nExp = 0,1,2,3,4,5 -> score = 0, 0.167, 0.333, 0.667, 0.833, 1.0
% This makes the threshold between 2 and 3 experimental neighbors explicit.
KTarget = max(1,round(KTarget));
bonusThreshold = min(max(round(bonusThreshold),1),KTarget);
nExp = min(max(round(nExp),0),KTarget);

if nExp >= bonusThreshold
    denom = max(KTarget - bonusThreshold + 1,1);
    frac = (nExp - bonusThreshold + 1) / denom;
    score = 0.5 + 0.5*frac.^powerVal;
else
    denom = max(bonusThreshold,1);
    frac = nExp / denom;
    score = 0.5*frac.^powerVal;
end

score = min(max(score,0),1);
end

function [nExp,coverageRatio,firstHasExp] = count_experimental_neighbors(expVec,ids,KTarget)
ids = ids(:)';
KTarget = max(1,round(KTarget));

nExp = 0;
coverageRatio = 0;
firstHasExp = false;

if isempty(ids) || isempty(expVec)
    return;
end

idsInt = round(ids);
validID = isfinite(idsInt) & idsInt >= 1 & idsInt <= numel(expVec);

if ~any(validID)
    return;
end

validExp = false(size(idsInt));
for k = 1:numel(idsInt)
    if validID(k)
        validExp(k) = is_not_missing(expVec(idsInt(k)));
    end
end

nExp = sum(validExp);
coverageRatio = min(nExp/KTarget,1);
firstHasExp = validExp(1);

end

function tf = is_not_missing(x)
% Robust missing-value test for numeric/cell/string experimental fields.

if iscell(x)
    if isempty(x) || isempty(x{1})
        tf = false;
        return;
    end
    x = x{1};
end

if isnumeric(x)
    tf = ~all(isnan(x(:))) && all(isfinite(x(~isnan(x))));
elseif isstring(x)
    s = strtrim(x);
    tf = strlength(s) > 0 && ~strcmpi(s,'NA') && ~strcmpi(s,'NaN');
elseif ischar(x)
    s = strtrim(x);
    tf = ~isempty(s) && ~strcmpi(s,'NA') && ~strcmpi(s,'NaN');
else
    tf = ~isempty(x);
end

end

function expVec = get_endpoint_exp_vector(CATMOS,endpoint)
switch endpoint
    case 'VT'
        expVec = CATMOS.model_VT.set.class_Exp;
    case 'NT'
        expVec = CATMOS.model_NT.set.class_Exp;
    case 'EPA'
        expVec = CATMOS.model_EPA.set.class_Exp;
    case 'GHS'
        expVec = CATMOS.model_GHS.set.class_Exp;
    case 'LD50'
        expVec = CATMOS.model_LD50.set.y_Exp_nAll;
    otherwise
        expVec = [];
end

if ~isempty(expVec) && size(expVec,2) > 1
    expVec = expVec(:,1);
end

end

function row = get_row(A,i)
if isempty(A)
    row = [];
elseif isvector(A)
    row = A(:)';
else
    row = A(i,:);
end
end

function val = get_res_scalar(res,fieldName,i,defaultVal)
if isfield(res,fieldName) && numel(res.(fieldName)) >= i
    val = res.(fieldName)(i,1);
else
    val = defaultVal;
end
end

function y = clip01(x)
if isempty(x) || ~isfinite(x)
    y = 0;
else
    y = min(max(x,0),1);
end
end

function res = store_diag_scalar(res,fieldName,i,value,templateField)
% Store scalar diagnostics without requiring external preallocation.
% If a template field exists, initialize diagnostic fields to matching length.
if ~isfield(res,fieldName)
    if isfield(res,templateField)
        n = size(res.(templateField),1);
    else
        n = i;
    end
    res.(fieldName) = NaN(n,1);
end

if size(res.(fieldName),1) < i
    res.(fieldName)(i,1) = NaN;
end
res.(fieldName)(i,1) = value;
end
