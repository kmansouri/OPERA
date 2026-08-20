function [res, diagOut] = catmos_conf_index_adjuster_v2(res,i,CATMOS,neighborsIDs,neighborsW,varargin)
%CATMOS_CONF_INDEX_ADJUSTER_V2 Post-WoE neighbor-based Conf_index adjustment.
%
% This helper is intended to be called AFTER woe_corr has selected and
% harmonized the CATMoS consensus prediction and AFTER consensus neighbors
% have been selected.
%
% Required inputs:
%   res         : CATMoS result structure containing Conf_index_CATMoS and
%                 CATMoS_LD50_pred for row i.
%   i           : row index.
%   CATMOS      : CATMOS object/structure containing CATMOS.model_LD50.set.
%   neighborsIDs: NxK matrix or 1xK vector of selected consensus neighbor IDs.
%   neighborsW  : NxK matrix or 1xK vector of selected neighbor raw weights.
%
% Common optional inputs:
%   'NeighborScore'       : NxK matrix or 1xK vector. If provided and
%                           WeightMode='auto' or 'score', this is used as
%                           the primary neighbor weight.
%   'WeightMode'          : 'auto' (default), 'score', 'rawW', or 'uniform'.
%   'UseCoverageWeight'   : true/false. Default true.
%   'CoverageMode'        : 'blend' (default), 'penalty', 'bonus', or 'none'.
%   'UsePredExpAgreement' : true/false. Default false, because this term only
%                           makes sense if YPredAll is truly predicted values
%                           on the same scale as y_Exp_nAll.
%   'YPredAll'            : vector of neighbor predicted/reference values.
%                           If empty and UsePredExpAgreement=true, the helper
%                           tries CATMOS.model_LD50.set.y.
%
% Main output:
%   res.Conf_index_CATMoS(i,1) is updated and clipped to [MinConf, MaxConf].
%
% Diagnostics stored in res by default:
%   Conf_index_CATMoS_preNN
%   Std_index_CATMoS
%   LocalY_index_CATMoS
%   PredExp_index_CATMoS
%   NN_exp_coverage_CATMoS
%   NN_valid_exp_count_CATMoS
%
% Notes:
%   All standard-deviation-derived terms are clipped to [0,1]. This prevents
%   1 - local_std/global_std from becoming negative when the selected
%   consensus neighbors are more spread than the global LD50 response.

p = inputParser;
p.addParameter('WeightMode','auto',@(x) ischar(x) || isstring(x));
p.addParameter('NeighborScore',[],@(x) isnumeric(x) || isempty(x));
p.addParameter('YPredAll',[],@(x) isnumeric(x) || isempty(x));

% Blend weights. If UsePredExpAgreement=false, PredExpWeight is set to 0 and
% the remaining weights are normalized when NormalizeBlendWeights=true.
p.addParameter('BaseWeight',0.50,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('ExpSpreadWeight',0.25,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('LocalYSpreadWeight',0.15,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('PredExpWeight',0.10,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('NormalizeBlendWeights',true,@(x) islogical(x) && isscalar(x));

% Coverage control based on how many of the selected neighbors have
% experimental LD50 values.
p.addParameter('UseCoverageWeight',true,@(x) islogical(x) && isscalar(x));
p.addParameter('CoverageMode','blend',@(x) ischar(x) || isstring(x));
p.addParameter('CoverageTarget',5,@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('CoveragePower',1,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('CoveragePenaltyStrength',0.25,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('CoverageBonusStrength',0.05,@(x) isnumeric(x) && isscalar(x) && x>=0);

% Term-specific behavior.
p.addParameter('UseOneExpNeighborCheck',true,@(x) islogical(x) && isscalar(x));
p.addParameter('UsePredExpAgreement',true,@(x) islogical(x) && isscalar(x));
p.addParameter('PredExpMetric','mae',@(x) ischar(x) || isstring(x));
p.addParameter('NeutralIfMissing',true,@(x) islogical(x) && isscalar(x));

% Final safety controls.
p.addParameter('PenaltyOnly',false,@(x) islogical(x) && isscalar(x));
p.addParameter('MaxBoost',0.10,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MaxDrop',0.30,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MinConf',0.10,@(x) isnumeric(x) && isscalar(x));
p.addParameter('MaxConf',1.00,@(x) isnumeric(x) && isscalar(x));
p.addParameter('StoreDiagnostics',false,@(x) islogical(x) && isscalar(x));

p.parse(varargin{:});
opt = p.Results;

baseConf = get_res_scalar(res,'Conf_index_CATMoS',i,NaN);
baseConf = clip01(baseConf);
if isnan(baseConf)
    baseConf = 0.5;
end

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

% Access LD50 experimental and reference/predicted vectors.
yExpAll = get_ld50_vector(CATMOS,'y_Exp_nAll');
yRefAll = get_ld50_vector(CATMOS,'y');

if ~isempty(opt.YPredAll)
    yPredAll = opt.YPredAll(:);
elseif opt.UsePredExpAgreement
    yPredAll = yRefAll;
else
    yPredAll = [];
end

% Validate IDs against yExp/yRef lengths.
idsInt = round(ids);
validNumericID = isfinite(idsInt) & idsInt >= 1;
validExpID = validNumericID & idsInt <= numel(yExpAll);
validRefID = validNumericID & idsInt <= numel(yRefAll);
validPredID = validNumericID & idsInt <= numel(yPredAll);

% Global scales.
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
% Term 1: experimental-neighbor spread index.
% Higher = selected neighbors with experimental LD50 values are internally
% more consistent. For one experimental neighbor, optionally compare the
% harmonized predicted LD50 against that neighbor's experimental value.
% -------------------------------------------------------------------------
yExp = NaN(1,K);
yExp(validExpID) = yExpAll(idsInt(validExpID));
validExp = validExpID & isfinite(yExp) & ~isnan(yExp);
nValidExp = sum(validExp);

if nValidExp >= 2
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
% Term 2: local y/reference spread index.
% This mirrors your old localY term based on CATMOS.model_LD50.set.y.
% Higher = the selected neighbors are locally consistent in y/reference space.
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
% Term 3: predicted-vs-experimental neighbor agreement.
% Use only if yPredAll is truly predicted values on the same scale as yExpAll.
% Higher = neighbors' predicted/reference values match their experimental values.
% -------------------------------------------------------------------------
PredExpIndex = neutral_value(baseConf,opt.NeutralIfMissing);
nValidPredExp = 0;

if opt.UsePredExpAgreement && ~isempty(yPredAll)
    yPred = NaN(1,K);
    yPred(validPredID) = yPredAll(idsInt(validPredID));
    validPredExp = validExp & validPredID & isfinite(yPred) & ~isnan(yPred);
    nValidPredExp = sum(validPredExp);

    if nValidPredExp >= 1
        e = yPred(validPredExp) - yExp(validPredExp);
        ww = normalize_weights(w(validPredExp));
        switch lower(char(opt.PredExpMetric))
            case 'rmse'
                err = sqrt(sum(ww .* (e.^2)));
            otherwise
                err = sum(ww .* abs(e));
        end
        PredExpIndex = 1 - err/globalStdExp;
    end
end
PredExpIndex = clip01(PredExpIndex);

% -------------------------------------------------------------------------
% Blend terms.
% -------------------------------------------------------------------------
baseW = opt.BaseWeight;
expW = opt.ExpSpreadWeight;
localW = opt.LocalYSpreadWeight;
predExpW = opt.PredExpWeight;

if ~opt.UsePredExpAgreement
    predExpW = 0;
end

if opt.NormalizeBlendWeights
    denom = baseW + expW + localW + predExpW;
    if denom > 0
        baseW = baseW/denom;
        expW = expW/denom;
        localW = localW/denom;
        predExpW = predExpW/denom;
    else
        baseW = 1; expW = 0; localW = 0; predExpW = 0;
    end
end

rawAdjustedConf = baseW*baseConf + ...
                  expW*ExpSpreadIndex + ...
                  localW*LocalYSpreadIndex + ...
                  predExpW*PredExpIndex;
rawAdjustedConf = clip01(rawAdjustedConf);

% Coverage factor based on how many selected neighbors have experimental LD50.
coverageFactor = min(nValidExp/opt.CoverageTarget,1)^opt.CoveragePower;
coverageFactor = clip01(coverageFactor);

if opt.UseCoverageWeight
    switch lower(char(opt.CoverageMode))
        case 'blend'
            candidateConf = baseConf + coverageFactor*(rawAdjustedConf - baseConf);
        case 'penalty'
            candidateConf = rawAdjustedConf * (1 - opt.CoveragePenaltyStrength*(1-coverageFactor));
        case 'bonus'
            candidateConf = rawAdjustedConf + opt.CoverageBonusStrength*coverageFactor;
        otherwise
            candidateConf = rawAdjustedConf;
    end
else
    candidateConf = rawAdjustedConf;
end

if opt.PenaltyOnly
    candidateConf = min(candidateConf,baseConf);
end

% Limit change relative to the WoE confidence to avoid extreme shifts.
candidateConf = min(candidateConf,baseConf + opt.MaxBoost);
candidateConf = max(candidateConf,baseConf - opt.MaxDrop);

finalConf = min(max(candidateConf,opt.MinConf),opt.MaxConf);

res.Conf_index_CATMoS(i,1) = finalConf;

if opt.StoreDiagnostics
    res.Conf_index_CATMoS_preNN(i,1) = baseConf;
    res.Std_index_CATMoS(i,1) = ExpSpreadIndex;
    res.LocalY_index_CATMoS(i,1) = LocalYSpreadIndex;
    res.PredExp_index_CATMoS(i,1) = PredExpIndex;
    res.NN_exp_coverage_CATMoS(i,1) = coverageFactor;
    res.NN_valid_exp_count_CATMoS(i,1) = nValidExp;
    res.NN_valid_predexp_count_CATMoS(i,1) = nValidPredExp;
end

diagOut = struct();
diagOut.baseConf = baseConf;
diagOut.finalConf = finalConf;
diagOut.rawAdjustedConf = rawAdjustedConf;
diagOut.ExpSpreadIndex = ExpSpreadIndex;
diagOut.LocalYSpreadIndex = LocalYSpreadIndex;
diagOut.PredExpIndex = PredExpIndex;
diagOut.nValidExp = nValidExp;
diagOut.nValidPredExp = nValidPredExp;
diagOut.coverageFactor = coverageFactor;
diagOut.weightsUsed = w;
diagOut.blendWeights = [baseW expW localW predExpW];
diagOut.globalStdExp = globalStdExp;
diagOut.globalStdRef = globalStdRef;

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
if sum(w) > 0
    w = w ./ sum(w);
else
    w = ones(size(w)) ./ numel(w);
end
end

function [mu,sigma] = weighted_mean_std(y,w)
y = y(:)';
w = normalize_weights(w);
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
