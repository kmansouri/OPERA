function [res, diagOut] = catmos_conf_index_adjuster_v10(res,i,CATMOS,neighborsIDs,neighborsW,varargin)
%CATMOS_CONF_INDEX_ADJUSTER_V10 Post-WoE neighbor-based Conf_index adjustment.
%
% Design requested in V8/V9/V10:
%   The final confidence is a weighted blend of:
%       BaseConf             = existing res.Conf_index_CATMoS from woe_corr
%       ExpSpreadIndex       = sparsity/spread of experimental LD50 neighbors
%       LocalYSpreadIndex    = sparsity/spread of neighbor predicted/reference y
%       PredExpIndex         = neighbor predicted/reference vs experimental agreement
%       CoverageEffectIndex  = availability-only coverage score
%
%   Default weights:
%       BaseWeight           = 0.50
%       ExpSpreadWeight      = 0.15
%       LocalYSpreadWeight   = 0.10
%       PredExpWeight        = 0.10
%       CoverageEffectWeight = 0.15
%
%   CoverageEffectIndex is high when >= CoverageThreshold neighbors have
%   experimental LD50 values, and low when < CoverageThreshold neighbors do.
%   With default CoverageThreshold=3 and CoverageTarget=5:
%       nExp = 0 -> 0.000
%       nExp = 1 -> 0.167
%       nExp = 2 -> 0.333
%       nExp = 3 -> 0.667
%       nExp = 4 -> 0.833
%       nExp = 5 -> 1.000
%
%   If nExp is 0 or 1, the ExpSpreadWeight is redistributed equally to
%   LocalYSpreadWeight, PredExpWeight, and CoverageEffectWeight by default.
%
%   V10 small edit:
%   If nExp is below CoverageThreshold, but the final CATMoS AD_index is high
%   and the FIRST selected consensus neighbor has experimental LD50, the
%   coverage component is made neutral instead of penalizing confidence.
%   This is controlled by:
%       UseHighADFirstExpNoCoveragePenalty = true
%       HighADFirstExpADIndexThreshold = 0.80
%
% This helper is intended to be called AFTER woe_corr and AFTER consensus
% nearest neighbors have been selected. It does not alter WoE.
%
% Recommended call:
%   [res,diag] = catmos_conf_index_adjuster_v10(res,i,CATMOS,neighborsIDs,neighborsW, ...
%       'NeighborScore',neighborsScore, ...
%       'WeightMode','auto');

p = inputParser;
p.addParameter('WeightMode','auto',@(x) ischar(x) || isstring(x));
p.addParameter('NeighborScore',[],@(x) isnumeric(x) || isempty(x));
p.addParameter('YPredAll',[],@(x) isnumeric(x) || isempty(x));

% Core component weights. These are normalized to sum to 1 by default after
% any redistribution.
p.addParameter('BaseWeight',0.50,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('ExpSpreadWeight',0.15,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('LocalYSpreadWeight',0.10,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('PredExpWeight',0.10,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('CoverageEffectWeight',0.15,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('NormalizeComponentWeights',true,@(x) islogical(x) && isscalar(x));

% Coverage score: availability-only score that is high for >= threshold and
% low for < threshold. This enters the blend as CoverageEffectIndex.
p.addParameter('CoverageTarget',5,@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('CoverageThreshold',3,@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('CoverageBonusPower',1.00,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('CoveragePenaltyPower',1.00,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('CoverageScoreMode','smooth_threshold',@(x) ischar(x) || isstring(x));

% V10 low-coverage exception:
% If coverage is low (< CoverageThreshold), but the final consensus AD_index
% is high and the FIRST selected consensus neighbor has experimental LD50,
% the availability component is made neutral instead of penalizing confidence.
p.addParameter('UseHighADFirstExpNoCoveragePenalty',true,@(x) islogical(x) && isscalar(x));
p.addParameter('HighADFirstExpADIndexThreshold',0.80,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);

% If very few experimental values are available, ExpSpreadIndex is weak or
% not meaningful. By default, when nExp is 0 or 1, move ExpSpreadWeight
% equally to LocalY, PredExp, and CoverageEffect.
p.addParameter('RedistributeExpWeightWhenSparse',true,@(x) islogical(x) && isscalar(x));
p.addParameter('ExpSparseRedistributeMaxN',1,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('RedistributeToLocalY',true,@(x) islogical(x) && isscalar(x));
p.addParameter('RedistributeToPredExp',true,@(x) islogical(x) && isscalar(x));
p.addParameter('RedistributeToCoverage',true,@(x) islogical(x) && isscalar(x));

% Term behavior.
p.addParameter('UseOneExpNeighborCheck',true,@(x) islogical(x) && isscalar(x));
p.addParameter('MinExpNeighborsForSpread',2,@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('UsePredExpAgreement',true,@(x) islogical(x) && isscalar(x));

% PredExp agreement modes:
%   'bin_overlap' (default): compare predicted vs experimental LD50 by WoE
%       bin agreement, using +/- ExpUncertaintyMargin around experimental LD50.
%   'value_error': old continuous-error style, normalized by globalStdExp.
p.addParameter('PredExpMode','bin_overlap',@(x) ischar(x) || isstring(x));
p.addParameter('PredExpMetric','mae',@(x) ischar(x) || isstring(x));
p.addParameter('ExpUncertaintyMargin',0.25,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('PredExpExactBinScore',1.00,@(x) isnumeric(x) && isscalar(x));
p.addParameter('PredExpOverlapScore',0.75,@(x) isnumeric(x) && isscalar(x));
p.addParameter('PredExpAdjacentScore',0.45,@(x) isnumeric(x) && isscalar(x));
p.addParameter('PredExpTwoAwayScore',0.20,@(x) isnumeric(x) && isscalar(x));
p.addParameter('PredExpFarScore',0.00,@(x) isnumeric(x) && isscalar(x));
p.addParameter('NeutralIfMissing',true,@(x) islogical(x) && isscalar(x));

% Safety controls. These are not forced low-coverage caps; they just prevent
% large jumps from any post-processing helper.
p.addParameter('PenaltyOnly',false,@(x) islogical(x) && isscalar(x));
p.addParameter('MaxBoost',0.10,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MaxDrop',0.35,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MinConf',0.10,@(x) isnumeric(x) && isscalar(x));
p.addParameter('MaxConf',1.00,@(x) isnumeric(x) && isscalar(x));
p.addParameter('StoreDiagnostics',false,@(x) islogical(x) && isscalar(x));

p.parse(varargin{:});
opt = p.Results;

% -------------------------------------------------------------------------
% Base confidence from WoE consensus.
% -------------------------------------------------------------------------
baseConf = get_res_scalar(res,'Conf_index_CATMoS',i,NaN);
baseConf = clip01(baseConf);
if isnan(baseConf)
    baseConf = 0.5;
end

% -------------------------------------------------------------------------
% Neighbor IDs and weights.
% -------------------------------------------------------------------------
ids = get_row_or_vector(neighborsIDs,i);
ids = ids(:)';
K = numel(ids);

rawW = get_row_or_vector(neighborsW,i);
rawW = rawW(:)';
if isempty(rawW) || numel(rawW) ~= K
    rawW = ones(1,K);
end

scoreW = [];
if ~isempty(opt.NeighborScore)
    scoreW = get_row_or_vector(opt.NeighborScore,i);
    scoreW = scoreW(:)';
    if numel(scoreW) ~= K
        scoreW = [];
    end
end

w = choose_neighbor_weights(rawW,scoreW,opt.WeightMode,K);

% -------------------------------------------------------------------------
% LD50 vectors.
%   yExpAll: experimental LD50 vector, usually log10 scale.
%   yRefAll: local predicted/reference y vector used by the old localY term.
%   yPredAll: optional predicted values for pred-vs-exp comparison.
% -------------------------------------------------------------------------
yExpAll = get_ld50_vector(CATMOS,'y_Exp_nAll');
yRefAll = get_ld50_vector(CATMOS,'y');

if ~isempty(opt.YPredAll)
    yPredAll = opt.YPredAll(:);
elseif opt.UsePredExpAgreement
    yPredAll = yRefAll;
else
    yPredAll = [];
end

idsInt = round(ids);
validNumericID = isfinite(idsInt) & idsInt >= 1;
validExpID  = validNumericID & idsInt <= numel(yExpAll);
validRefID  = validNumericID & idsInt <= numel(yRefAll);
validPredID = validNumericID & idsInt <= numel(yPredAll);

% Global scales used to normalize spread/error indices.
globalStdExp = safe_std(yExpAll(isfinite(yExpAll) & ~isnan(yExpAll)));
if ~isfinite(globalStdExp) || globalStdExp <= 0
    globalStdExp = safe_std(yRefAll(isfinite(yRefAll) & ~isnan(yRefAll)));
end
if ~isfinite(globalStdExp) || globalStdExp <= 0
    globalStdExp = 1;
end

globalStdRef = safe_std(yRefAll(isfinite(yRefAll) & ~isnan(yRefAll)));
if ~isfinite(globalStdRef) || globalStdRef <= 0
    globalStdRef = globalStdExp;
end
if ~isfinite(globalStdRef) || globalStdRef <= 0
    globalStdRef = 1;
end

% -------------------------------------------------------------------------
% Experimental availability.
% -------------------------------------------------------------------------
yExp = NaN(1,K);
yExp(validExpID) = yExpAll(idsInt(validExpID));
validExp = validExpID & isfinite(yExp) & ~isnan(yExp);
nValidExp = sum(validExp);

coverageTarget = max(1,round(opt.CoverageTarget));
coverageThreshold = min(max(1,round(opt.CoverageThreshold)),coverageTarget);
coverageRatio = clip01(nValidExp/coverageTarget);
coverageEffectIndex = coverage_score(nValidExp,coverageTarget,coverageThreshold, ...
    opt.CoverageScoreMode,opt.CoverageBonusPower,opt.CoveragePenaltyPower);
coverageEffectIndexRaw = coverageEffectIndex;

% V10: low-coverage penalty exception for high-AD predictions with an
% experimental value on the first selected consensus neighbor. This affects
% only the availability/coverage component; ExpSpreadIndex and PredExpIndex
% still reflect the actual experimental evidence quality.
ADindexForCoverageGuard = get_res_scalar(res,'AD_index_CATMoS',i,NaN);
firstNNHasExp = false;
if K >= 1
    firstNNHasExp = validExp(1);
end
coveragePenaltyBypassed = false;
if opt.UseHighADFirstExpNoCoveragePenalty && ...
        nValidExp < coverageThreshold && ...
        firstNNHasExp && ...
        isfinite(ADindexForCoverageGuard) && ...
        ADindexForCoverageGuard >= opt.HighADFirstExpADIndexThreshold

    % Neutralize only the low-coverage penalty by setting the coverage score
    % to at least baseConf. If baseConf is high, the coverage component no
    % longer pulls confidence down; if baseConf is low, this does not create
    % an artificial high-coverage bonus.
    coverageEffectIndex = max(coverageEffectIndex,baseConf);
    coverageEffectIndex = clip01(coverageEffectIndex);
    coveragePenaltyBypassed = true;
end

% -------------------------------------------------------------------------
% Term A: ExpSpreadIndex. Experimental-dependent.
% Higher = selected neighbors with experimental LD50 are internally consistent.
% -------------------------------------------------------------------------
if nValidExp >= opt.MinExpNeighborsForSpread
    y = yExp(validExp);
    ww = normalize_weights(w(validExp));
    [~,localStdExp] = weighted_mean_std(y,ww);
    ExpSpreadIndex = 1 - localStdExp/globalStdExp;
elseif nValidExp == 1 && opt.UseOneExpNeighborCheck
    predLD50 = get_res_scalar(res,'CATMoS_LD50_pred',i,NaN);
    yOne = yExp(validExp);
    if isfinite(predLD50) && isfinite(yOne)
        ExpSpreadIndex = 1 - abs(predLD50 - yOne)/globalStdExp;
    else
        ExpSpreadIndex = neutral_value(baseConf,opt.NeutralIfMissing);
    end
else
    ExpSpreadIndex = neutral_value(baseConf,opt.NeutralIfMissing);
end
ExpSpreadIndex = clip01(ExpSpreadIndex);

% -------------------------------------------------------------------------
% Term B: LocalYSpreadIndex. Not experimental-dependent.
% Higher = selected neighbors are locally consistent in y/reference values.
% -------------------------------------------------------------------------
yRef = NaN(1,K);
yRef(validRefID) = yRefAll(idsInt(validRefID));
validRef = validRefID & isfinite(yRef) & ~isnan(yRef);

if sum(validRef) >= 2
    y = yRef(validRef);
    ww = normalize_weights(w(validRef));
    [~,localStdRef] = weighted_mean_std(y,ww);
    LocalYSpreadIndex = 1 - localStdRef/globalStdRef;
else
    LocalYSpreadIndex = neutral_value(baseConf,opt.NeutralIfMissing);
end
LocalYSpreadIndex = clip01(LocalYSpreadIndex);

% -------------------------------------------------------------------------
% Term C: PredExpIndex. Experimental-dependent.
% Higher = selected neighbors' predicted/reference values agree with exp LD50.
%
% V9 default uses a bin-overlap score instead of a raw continuous error:
%   1.00 = predicted and experimental central values are in the same WoE bin
%   0.75 = predicted bin overlaps the experimental +/- margin bin range
%   0.45 = predicted bin is one bin away from that experimental range
%   0.20 = predicted bin is two bins away
%   0.00 = farther away
% This is more aligned with Conf_index_CATMoS as confidence in the harmonized
% ordinal toxicity bin rather than exact continuous LD50 equality.
% -------------------------------------------------------------------------
PredExpIndex = neutral_value(baseConf,opt.NeutralIfMissing);
nValidPredExp = 0;

if opt.UsePredExpAgreement && ~isempty(yPredAll)
    yPred = NaN(1,K);
    yPred(validPredID) = yPredAll(idsInt(validPredID));
    validPredExp = validExp & validPredID & isfinite(yPred) & ~isnan(yPred);
    nValidPredExp = sum(validPredExp);

    if nValidPredExp >= 1
        yPredPE = yPred(validPredExp);
        yExpPE  = yExp(validPredExp);
        ww = normalize_weights(w(validPredExp));

        switch lower(char(opt.PredExpMode))
            case {'bin','bins','bin_overlap','binoverlap'}
                predExpScores = predexp_bin_overlap_scores(yPredPE,yExpPE, ...
                    'ExpMargin',opt.ExpUncertaintyMargin, ...
                    'ExactBinScore',opt.PredExpExactBinScore, ...
                    'OverlapScore',opt.PredExpOverlapScore, ...
                    'AdjacentScore',opt.PredExpAdjacentScore, ...
                    'TwoAwayScore',opt.PredExpTwoAwayScore, ...
                    'FarScore',opt.PredExpFarScore);
                PredExpIndex = sum(ww .* predExpScores);

            otherwise
                % Optional old-style continuous error mode. This remains
                % available for comparison/testing, but is no longer default.
                e = yPredPE - yExpPE;
                switch lower(char(opt.PredExpMetric))
                    case 'rmse'
                        err = sqrt(sum(ww .* (e.^2)));
                    otherwise
                        err = sum(ww .* abs(e));
                end
                PredExpIndex = 1 - err/globalStdExp;
        end
    end
end
PredExpIndex = clip01(PredExpIndex);

% -------------------------------------------------------------------------
% Component weights, including requested redistribution.
% -------------------------------------------------------------------------
baseW = opt.BaseWeight;
expW = opt.ExpSpreadWeight;
localW = opt.LocalYSpreadWeight;
predExpW = opt.PredExpWeight;
coverageW = opt.CoverageEffectWeight;

if ~opt.UsePredExpAgreement
    % PredExpIndex is unavailable by choice. Move its weight equally to
    % LocalY and CoverageEffect by default, then set PredExp weight to 0.
    moveW = predExpW;
    predExpW = 0;
    localW = localW + moveW/2;
    coverageW = coverageW + moveW/2;
end

if opt.RedistributeExpWeightWhenSparse && nValidExp <= opt.ExpSparseRedistributeMaxN && expW > 0
    moveW = expW;
    expW = 0;

    receivers = [];
    if opt.RedistributeToLocalY
        receivers(end+1) = 1; %#ok<AGROW>
    end
    if opt.RedistributeToPredExp && opt.UsePredExpAgreement
        receivers(end+1) = 2; %#ok<AGROW>
    end
    if opt.RedistributeToCoverage
        receivers(end+1) = 3; %#ok<AGROW>
    end

    if isempty(receivers)
        localW = localW + moveW;
    else
        addW = moveW/numel(receivers);
        for rr = receivers
            switch rr
                case 1
                    localW = localW + addW;
                case 2
                    predExpW = predExpW + addW;
                case 3
                    coverageW = coverageW + addW;
            end
        end
    end
end

if opt.NormalizeComponentWeights
    denom = baseW + expW + localW + predExpW + coverageW;
    if denom > 0
        baseW = baseW/denom;
        expW = expW/denom;
        localW = localW/denom;
        predExpW = predExpW/denom;
        coverageW = coverageW/denom;
    else
        baseW = 1; expW = 0; localW = 0; predExpW = 0; coverageW = 0;
    end
end

candidateConf = baseW*baseConf + ...
                expW*ExpSpreadIndex + ...
                localW*LocalYSpreadIndex + ...
                predExpW*PredExpIndex + ...
                coverageW*coverageEffectIndex;

candidateConf = clip01(candidateConf);

if opt.PenaltyOnly
    candidateConf = min(candidateConf,baseConf);
end

candidateConfBeforeCaps = candidateConf;

% Limit changes relative to the original WoE confidence. These are general
% safety controls, not low-coverage caps.
candidateConf = min(candidateConf,baseConf + opt.MaxBoost);
candidateConf = max(candidateConf,baseConf - opt.MaxDrop);

finalConf = min(max(candidateConf,opt.MinConf),opt.MaxConf);
res.Conf_index_CATMoS(i,1) = finalConf;

% -------------------------------------------------------------------------
% Diagnostics.
% -------------------------------------------------------------------------
if opt.StoreDiagnostics
    res.Conf_index_CATMoS_preNN(i,1) = baseConf;
    res.Std_index_CATMoS(i,1) = ExpSpreadIndex;
    res.LocalY_index_CATMoS(i,1) = LocalYSpreadIndex;
    res.PredExp_index_CATMoS(i,1) = PredExpIndex;
    res.NN_exp_coverage_ratio_CATMoS(i,1) = coverageRatio;
    res.NN_exp_coverage_score_CATMoS(i,1) = coverageEffectIndex;
    res.NN_exp_coverage_score_raw_CATMoS(i,1) = coverageEffectIndexRaw;
    res.NN_first_neighbor_has_exp_CATMoS(i,1) = double(firstNNHasExp);
    res.NN_AD_index_for_coverage_guard_CATMoS(i,1) = ADindexForCoverageGuard;
    res.NN_highAD_firstExp_no_coverage_penalty_CATMoS(i,1) = double(coveragePenaltyBypassed);
    res.NN_valid_exp_count_CATMoS(i,1) = nValidExp;
    res.NN_valid_predexp_count_CATMoS(i,1) = nValidPredExp;
    res.NN_conf_base_weight_CATMoS(i,1) = baseW;
    res.NN_conf_exp_weight_CATMoS(i,1) = expW;
    res.NN_conf_localY_weight_CATMoS(i,1) = localW;
    res.NN_conf_predexp_weight_CATMoS(i,1) = predExpW;
    res.NN_conf_coverage_weight_CATMoS(i,1) = coverageW;
    res.NN_conf_candidate_before_caps_CATMoS(i,1) = candidateConfBeforeCaps;
end

diagOut = struct();
diagOut.baseConf = baseConf;
diagOut.finalConf = finalConf;
diagOut.candidateConfBeforeCaps = candidateConfBeforeCaps;
diagOut.ExpSpreadIndex = ExpSpreadIndex;
diagOut.LocalYSpreadIndex = LocalYSpreadIndex;
diagOut.PredExpIndex = PredExpIndex;
diagOut.CoverageEffectIndex = coverageEffectIndex;
diagOut.CoverageEffectIndexRaw = coverageEffectIndexRaw;
diagOut.firstNNHasExp = firstNNHasExp;
diagOut.ADindexForCoverageGuard = ADindexForCoverageGuard;
diagOut.coveragePenaltyBypassed = coveragePenaltyBypassed;
diagOut.HighADFirstExpADIndexThreshold = opt.HighADFirstExpADIndexThreshold;
diagOut.nValidExp = nValidExp;
diagOut.nValidPredExp = nValidPredExp;
diagOut.coverageRatio = coverageRatio;
diagOut.weightsUsed = w;
diagOut.blendWeights = struct('baseW',baseW,'expW',expW,'localW',localW, ...
    'predExpW',predExpW,'coverageW',coverageW);
diagOut.globalStdExp = globalStdExp;
diagOut.globalStdRef = globalStdRef;
diagOut.PredExpMode = char(opt.PredExpMode);
diagOut.ExpUncertaintyMargin = opt.ExpUncertaintyMargin;

end

% =============================== helpers =================================

function v = get_row_or_vector(A,i)
if isempty(A)
    v = [];
elseif isvector(A)
    v = A;
else
    v = A(i,:);
end
end

function val = get_res_scalar(res,field,i,defaultVal)
if isfield(res,field) && size(res.(field),1) >= i && ~isempty(res.(field)(i,1))
    val = res.(field)(i,1);
else
    val = defaultVal;
end
end

function y = get_ld50_vector(CATMOS,field)
y = [];
try
    if isfield(CATMOS,'model_LD50') && isfield(CATMOS.model_LD50,'set') && isfield(CATMOS.model_LD50.set,field)
        y = CATMOS.model_LD50.set.(field)(:);
    end
catch
    y = [];
end
if isempty(y)
    y = NaN(0,1);
end
end

function w = choose_neighbor_weights(rawW,scoreW,mode,K)
mode = lower(char(mode));
switch mode
    case 'score'
        if ~isempty(scoreW)
            w = scoreW;
        else
            w = rawW;
        end
    case 'raww'
        w = rawW;
    case 'uniform'
        w = ones(1,K);
    otherwise % auto
        if ~isempty(scoreW) && any(isfinite(scoreW) & scoreW>0)
            w = scoreW;
        else
            w = rawW;
        end
end
w = w(:)';
if isempty(w) || numel(w) ~= K
    w = ones(1,K);
end
w(~isfinite(w)) = 0;
w = max(w,0);
if sum(w) <= 0
    w = ones(1,K);
end
end

function w = normalize_weights(w)
w = w(:)';
w(~isfinite(w)) = 0;
w = max(w,0);
if isempty(w)
    w = [];
elseif sum(w) > 0
    w = w ./ sum(w);
else
    w = ones(size(w)) ./ numel(w);
end
end

function [mu,sigma] = weighted_mean_std(y,w)
y = y(:)';
w = normalize_weights(w);
if isempty(y) || isempty(w)
    mu = NaN;
    sigma = NaN;
    return;
end
mu = sum(w .* y);
sigma = sqrt(sum(w .* (y-mu).^2));
end

function s = safe_std(x)
x = x(:);
x = x(isfinite(x) & ~isnan(x));
if numel(x) >= 2
    s = std(x);
else
    s = NaN;
end
end

function y = clip01(x)
y = x;
y(~isfinite(y)) = NaN;
y = min(max(y,0),1);
end

function v = neutral_value(baseConf,useBase)
if useBase
    v = baseConf;
else
    v = 0.5;
end
end


function scores = predexp_bin_overlap_scores(yPred,yExp,varargin)
%PREDEXP_BIN_OVERLAP_SCORES Bin-based PredExp agreement for log10 LD50.
%
% yPred and yExp must be log10(LD50) values.
% Score rules:
%   same central bin                        -> ExactBinScore
%   predicted bin overlaps yExp +/- margin  -> OverlapScore
%   one bin away from uncertainty range      -> AdjacentScore
%   two bins away                            -> TwoAwayScore
%   farther                                  -> FarScore
%
% CATMoS/WoE LD50 bins:
%   1 <=5, 2 >5-50, 3 >50-300, 4 >300-500,
%   5 >500-2000, 6 >2000-5000, 7 >5000 mg/kg.

p = inputParser;
p.addParameter('ExpMargin',0.25,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('ExactBinScore',1.00,@(x) isnumeric(x) && isscalar(x));
p.addParameter('OverlapScore',0.75,@(x) isnumeric(x) && isscalar(x));
p.addParameter('AdjacentScore',0.45,@(x) isnumeric(x) && isscalar(x));
p.addParameter('TwoAwayScore',0.20,@(x) isnumeric(x) && isscalar(x));
p.addParameter('FarScore',0.00,@(x) isnumeric(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

yPred = yPred(:)';
yExp  = yExp(:)';
n = min(numel(yPred),numel(yExp));
yPred = yPred(1:n);
yExp = yExp(1:n);

edges = [-Inf, log10([5,50,300,500,2000,5000]), Inf];
scores = zeros(1,n);

for k = 1:n
    if ~isfinite(yPred(k)) || ~isfinite(yExp(k))
        scores(k) = NaN;
        continue;
    end

    predBin = ld50_log_to_woe_bin(yPred(k),edges);
    expBin  = ld50_log_to_woe_bin(yExp(k),edges);

    expLo = yExp(k) - opt.ExpMargin;
    expHi = yExp(k) + opt.ExpMargin;

    % Bins whose ranges overlap the experimental uncertainty interval.
    expUncBins = find(edges(1:end-1) < expHi & edges(2:end) >= expLo);
    if isempty(expUncBins)
        expUncBins = expBin;
    end

    if predBin == expBin
        scores(k) = opt.ExactBinScore;
    elseif ismember(predBin,expUncBins)
        scores(k) = opt.OverlapScore;
    else
        gap = min(abs(predBin - expUncBins));
        if gap == 1
            scores(k) = opt.AdjacentScore;
        elseif gap == 2
            scores(k) = opt.TwoAwayScore;
        else
            scores(k) = opt.FarScore;
        end
    end
end

scores(~isfinite(scores)) = 0;
scores = clip01(scores);
end

function bin = ld50_log_to_woe_bin(x,edges)
%LD50_LOG_TO_WOE_BIN Convert log10 LD50 to WoE bin 1..7.
% Uses <= upper edge so exact thresholds stay in the lower/more toxic bin.
bin = find(x <= edges(2:end),1,'first');
if isempty(bin)
    bin = 7;
end
end

function c = coverage_score(nValidExp,coverageTarget,coverageThreshold,mode,bonusPower,penaltyPower)
mode = lower(char(mode));
n = min(max(nValidExp,0),coverageTarget);
T = coverageThreshold;
K = coverageTarget;

switch mode
    case 'binary_threshold'
        if n >= T
            c = 1;
        else
            c = 0;
        end
    case 'linear_ratio'
        c = n/K;
    otherwise % smooth_threshold
        if n >= T
            % n=T starts above neutral; n=K reaches 1.
            denom = max(K - T + 1,eps);
            bonusFrac = (n - T + 1)/denom;
            c = 0.5 + 0.5*(bonusFrac^bonusPower);
        else
            % n<T stays below neutral. For T=3: n=0,1,2 -> 0,.167,.333.
            penaltyFrac = n/T;
            c = 0.5*(penaltyFrac^penaltyPower);
        end
end
c = clip01(c);
end
