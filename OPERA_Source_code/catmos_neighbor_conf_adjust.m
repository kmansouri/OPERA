function [res, Std_index_CATMoS, localY_index_CATMoS, nnConfDiag] = catmos_neighbor_conf_adjust(res,i,CATMOS,neighborsIDs,neighborsW,varargin)
%CATMOS_NEIGHBOR_STD_CONF_ADJUST Post-WoE neighbor-STD confidence adjustment.
%
% This helper reproduces the intent of the original CATMoS neighbor-based
% confidence adjustment, but makes the STD-derived indices bounded in [0,1]
% and optionally weights the adjustment by the number of selected neighbors
% with experimental LD50 values.
%
% Typical use after woe_corr and after selecting consensus neighbors:
%
%   [res, Std_index_CATMoS, localY_index_CATMoS, diag] = ...
%       catmos_neighbor_std_conf_adjust(res,i,CATMOS,neighborsIDs,neighborsW);
%
% Inputs:
%   res          : result structure containing res.Conf_index_CATMoS and
%                  res.CATMoS_LD50_pred for row i.
%   i            : row index.
%   CATMOS       : CATMOS model structure. Uses:
%                    CATMOS.model_LD50.set.y_Exp_nAll
%                    CATMOS.model_LD50.set.y
%   neighborsIDs : N x K matrix or 1 x K vector of selected neighbor IDs.
%   neighborsW   : N x K matrix or 1 x K vector of selected neighbor weights.
%
% Outputs:
%   res                      : updated res with adjusted Conf_index_CATMoS.
%   Std_index_CATMoS          : bounded [0,1] experimental-neighbor STD index.
%   localY_index_CATMoS       : bounded [0,1] training-y local STD index.
%   nnConfDiag                : diagnostic structure.
%
% Key options:
%   'UseCoverageWeight'       : true/false. Default true.
%   'CoverageMode'            : 'blend', 'penalty', 'bonus', or 'none'.
%   'ExpectedNeighborCount'   : default 5.
%   'CoveragePower'           : default 1. Higher makes low coverage matter more.
%   'BaseWeight'              : default 0.5.
%   'StdWeight'               : default 0.3.
%   'LocalYWeight'            : default 0.2.
%   'MinConf'                 : default 0.1.
%
% Recommended default behaviour:
%   UseCoverageWeight=true and CoverageMode='blend'. This means when few or no
%   selected neighbors have experimental LD50 values, the helper applies less
%   of the neighbor-based adjustment and keeps the original WoE confidence.

p = inputParser;
p.addParameter('UseCoverageWeight',true,@(x) islogical(x) && isscalar(x));
p.addParameter('CoverageMode','blend',@(x) ischar(x) || isstring(x));
p.addParameter('ExpectedNeighborCount',5,@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('CoveragePower',1,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('CoveragePenaltyStrength',0.25,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('CoverageBonusStrength',0.05,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('BaseWeight',0.5,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('StdWeight',0.3,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('LocalYWeight',0.2,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MinConf',0.1,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('MaxConf',1.0,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('NormalizeBlendWeights',false,@(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
coverageMode = lower(char(opt.CoverageMode));

% -------------------------------------------------------------------------
% Pull the neighbor row and clean basic inputs.
% -------------------------------------------------------------------------
ids = get_row_or_vector(neighborsIDs,i);
w   = get_row_or_vector(neighborsW,i);

ids = ids(:)';
w   = w(:)';

K = min(numel(ids),numel(w));
ids = ids(1:K);
w = w(1:K);

baseConf = get_res_scalar(res,'Conf_index_CATMoS',i,NaN);
baseConf = clip01(baseConf);

% If the required CATMOS LD50 fields are missing, keep base confidence.
if ~isfield(CATMOS,'model_LD50') || ~isfield(CATMOS.model_LD50,'set') || ...
        ~isfield(CATMOS.model_LD50.set,'y_Exp_nAll') || ...
        ~isfield(CATMOS.model_LD50.set,'y')
    Std_index_CATMoS = NaN;
    localY_index_CATMoS = NaN;
    nnConfDiag = empty_diag(baseConf,NaN,NaN,0,0,NaN,baseConf,'missing CATMOS LD50 fields');
    res.Conf_index_CATMoS(i,1) = baseConf;
    return;
end

yExpAll = CATMOS.model_LD50.set.y_Exp_nAll(:);
yTrain  = CATMOS.model_LD50.set.y(:);

validID_exp = isfinite(ids) & ids>=1 & ids<=numel(yExpAll) & abs(ids-round(ids))<1e-9;
validID_y   = isfinite(ids) & ids>=1 & ids<=numel(yTrain)  & abs(ids-round(ids))<1e-9;

idxExp = NaN(size(ids));
idxY   = NaN(size(ids));
idxExp(validID_exp) = round(ids(validID_exp));
idxY(validID_y)     = round(ids(validID_y));

% Experimental LD50 values for selected neighbors.
yExp = NaN(size(ids));
if any(validID_exp)
    yExp(validID_exp) = yExpAll(idxExp(validID_exp));
end
validExp = validID_exp & isfinite(yExp) & ~isnan(yExp);

% Training/model y values for selected neighbors.
yLocal = NaN(size(ids));
if any(validID_y)
    yLocal(validID_y) = yTrain(idxY(validID_y));
end
validY = validID_y & isfinite(yLocal) & ~isnan(yLocal) & isfinite(w);

nValidExp = sum(validExp);
expectedK = opt.ExpectedNeighborCount;
coverageRaw = min(nValidExp/expectedK,1);
coverageFactor = coverageRaw.^opt.CoveragePower;

% Global standard deviations.
globalStdExp = std(yExpAll(isfinite(yExpAll) & ~isnan(yExpAll)));
globalStdY   = std(yTrain(isfinite(yTrain) & ~isnan(yTrain)));

% -------------------------------------------------------------------------
% Term 1: Std_index_CATMoS, based on experimental LD50 values of neighbors.
% This follows the structure of the original code, but bounds the output.
% -------------------------------------------------------------------------
Std_index_CATMoS = 0;

firstIDValidExp = false;
firstYExp = NaN;
firstW = NaN;
if K >= 1 && validExp(1)
    firstIDValidExp = true;
    firstYExp = yExp(1);
    firstW = w(1);
end

if ~isfinite(globalStdExp) || globalStdExp <= 0
    Std_index_CATMoS = baseConf;

elseif K >= 1 && isfinite(w(1)) && w(1)==1 && firstIDValidExp
    Std_index_CATMoS = 1;

elseif ~firstIDValidExp && nValidExp <= 1
    Std_index_CATMoS = 0;

elseif firstIDValidExp && nValidExp <= 1
    predLD50 = get_res_scalar(res,'CATMoS_LD50_pred',i,NaN);
    if isfinite(predLD50)
        localStd = local_weighted_std([predLD50, firstYExp],[1, max(firstW,0)]);
        Std_index_CATMoS = 1 - localStd/globalStdExp;
    else
        Std_index_CATMoS = 0;
    end

else
    yExpValid = yExp(validExp);
    wExpValid = w(validExp);
    wExpValid(~isfinite(wExpValid)) = 0;
    wExpValid = max(wExpValid,0);

    if sum(wExpValid) == 0
        localStd = std(yExpValid);
    else
        localStd = local_weighted_std(yExpValid,wExpValid);
    end
    Std_index_CATMoS = 1 - localStd/globalStdExp;
end

Std_index_CATMoS = clip01(Std_index_CATMoS);

% -------------------------------------------------------------------------
% Term 2: localY_index_CATMoS, based on local spread in CATMOS.model_LD50.set.y.
% -------------------------------------------------------------------------
if ~isfinite(globalStdY) || globalStdY <= 0
    localY_index_CATMoS = baseConf;

elseif sum(validY) >= 2
    yValid = yLocal(validY);
    wValid = w(validY);
    wValid(~isfinite(wValid)) = 0;
    wValid = max(wValid,0);

    if sum(wValid) == 0
        localStdY = std(yValid);
    else
        localStdY = local_weighted_std(yValid,wValid);
    end
    localY_index_CATMoS = 1 - localStdY/globalStdY;

elseif sum(validY) == 1
    % One available neighbor has no local spread, so keep this term high.
    localY_index_CATMoS = 1;

else
    % No usable y values: do not create an artificial penalty.
    localY_index_CATMoS = baseConf;
end

localY_index_CATMoS = clip01(localY_index_CATMoS);

% -------------------------------------------------------------------------
% Raw old-style blend, with bounded components.
% -------------------------------------------------------------------------
baseW  = opt.BaseWeight;
stdW   = opt.StdWeight;
localW = opt.LocalYWeight;

if opt.NormalizeBlendWeights
    denom = baseW + stdW + localW;
    if denom > 0
        baseW = baseW/denom;
        stdW = stdW/denom;
        localW = localW/denom;
    end
end

rawAdjustedConf = baseW*baseConf + stdW*Std_index_CATMoS + localW*localY_index_CATMoS;
rawAdjustedConf = clip01(rawAdjustedConf);

% -------------------------------------------------------------------------
% Optional coverage weighting.
%   blend  : coverage controls how much the neighbor adjustment moves baseConf.
%   penalty: gradually penalize adjusted confidence if fewer than expectedK
%            experimental neighbors are available.
%   bonus  : gradually boost adjusted confidence when coverage is high.
%   none   : do not use coverage.
% -------------------------------------------------------------------------
if opt.UseCoverageWeight
    switch coverageMode
        case 'blend'
            finalConf = baseConf + coverageFactor*(rawAdjustedConf - baseConf);

        case 'penalty'
            coveragePenalty = 1 - opt.CoveragePenaltyStrength*(1-coverageFactor);
            finalConf = rawAdjustedConf * coveragePenalty;

        case 'bonus'
            finalConf = rawAdjustedConf + opt.CoverageBonusStrength*coverageFactor;

        case 'none'
            finalConf = rawAdjustedConf;

        otherwise
            error('Unknown CoverageMode: %s. Use blend, penalty, bonus, or none.',coverageMode);
    end
else
    finalConf = rawAdjustedConf;
end

finalConf = min(max(finalConf,opt.MinConf),opt.MaxConf);
res.Conf_index_CATMoS(i,1) = finalConf;

% Optional diagnostics. Store in res if these fields are useful downstream.
%res.Std_index_CATMoS(i,1) = Std_index_CATMoS;
%res.LocalY_index_CATMoS(i,1) = localY_index_CATMoS;
%res.NN_exp_coverage_CATMoS(i,1) = coverageFactor;
%res.NN_valid_exp_count_CATMoS(i,1) = nValidExp;

nnConfDiag = struct();
nnConfDiag.baseConf = baseConf;
nnConfDiag.Std_index_CATMoS = Std_index_CATMoS;
nnConfDiag.localY_index_CATMoS = localY_index_CATMoS;
nnConfDiag.rawAdjustedConf = rawAdjustedConf;
nnConfDiag.finalConf = finalConf;
nnConfDiag.nValidExp = nValidExp;
nnConfDiag.expectedNeighborCount = expectedK;
nnConfDiag.coverageRaw = coverageRaw;
nnConfDiag.coverageFactor = coverageFactor;
nnConfDiag.coverageMode = coverageMode;
nnConfDiag.globalStdExp = globalStdExp;
nnConfDiag.globalStdY = globalStdY;
nnConfDiag.validExpMask = validExp;
nnConfDiag.validYMask = validY;

end

%% ========================================================================
% Helper functions
%% ========================================================================
function row = get_row_or_vector(X,i)
if isempty(X)
    row = [];
elseif isvector(X)
    row = X(:)';
else
    if size(X,1) >= i
        row = X(i,:);
    else
        row = [];
    end
end
end

function val = get_res_scalar(res,fieldName,i,defaultVal)
if isfield(res,fieldName) && numel(res.(fieldName)) >= i
    val = res.(fieldName)(i,1);
else
    val = defaultVal;
end
end

function x = clip01(x)
x(~isfinite(x)) = 0;
x = min(max(x,0),1);
end

function s = local_weighted_std(x,w)
x = x(:);
w = w(:);
valid = isfinite(x) & isfinite(w) & w>=0;
x = x(valid);
w = w(valid);

if isempty(x)
    s = NaN;
    return;
end

if numel(x) == 1
    s = 0;
    return;
end

if sum(w) <= 0
    s = std(x);
    return;
end

w = w ./ sum(w);
mu = sum(w.*x);
s = sqrt(sum(w.*(x-mu).^2));
end

function d = empty_diag(baseConf,Std_index,localY_index,nValidExp,expectedK,coverageFactor,finalConf,msg)
d = struct();
d.baseConf = baseConf;
d.Std_index_CATMoS = Std_index;
d.localY_index_CATMoS = localY_index;
d.rawAdjustedConf = baseConf;
d.finalConf = finalConf;
d.nValidExp = nValidExp;
d.expectedNeighborCount = expectedK;
d.coverageRaw = NaN;
d.coverageFactor = coverageFactor;
d.coverageMode = '';
d.globalStdExp = NaN;
d.globalStdY = NaN;
d.validExpMask = [];
d.validYMask = [];
d.message = msg;
end
