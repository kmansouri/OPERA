function [neighborsIDs, neighborsSrc, neighborsScore, neighborsDC, neighborsW] = catmos_consensus_neighbors_with_weights(res,i,predVT,predNT,predEPA,predGHS,predLD50,varargin)
%CATMOS_CONSENSUS_NEIGHBORS_WITH_WEIGHTS Select final-consensus neighbors.
%
% Usage:
%   [neighborsIDs,neighborsSrc,neighborsScore,neighborsDC,neighborsW] = ...
%       catmos_consensus_neighbors_with_weights(res,i,predVT,predNT,predEPA,predGHS,predLD50);
%
% Required:
%   res.CATMoS_winning_model_names{i,1} must be a cell array such as
%   {'NT','GHS'} or {'EPA','GHS','LD50'}.
%
% Outputs:
%   neighborsIDs    : 1 x K final selected neighbor IDs
%   neighborsSrc    : 1 x K cell array showing source model(s), e.g. 'GHS' or 'EPA+LD50'
%   neighborsScore  : 1 x K final cross-endpoint consensus ranking score
%   neighborsDC     : 1 x K source-model Euclidean distance for the selected occurrence
%   neighborsW      : 1 x K raw source-model nearest-neighbor weight for the selected occurrence
%
% Notes:
%   neighborsW is the raw predX.w value associated with the selected neighbor.
%   If the final selected neighbors all come from LD50 and no duplicate merging changes
%   anything, neighborsW is predLD50.w(i,:) reordered to match neighborsIDs.
%   If a neighbor appears from multiple endpoint models, this function keeps the
%   occurrence with the highest final candidate score and stores that occurrence's raw weight.
%
% Options:
%   'K'                         default 5
%   'Balanced'                  default true; select one best neighbor per winning model first
%   'UseWeights'                default true; use predX.w if available
%   'WeightComponent'           default 0.60
%   'LocalDistanceComponent'    default 0.25
%   'GlobalDistanceComponent'   default 0.15
%   'ADGlobalFloor'             default 0.25
%   'ReliabilityFloor'          default 0.05

p = inputParser;
p.addParameter('K',5,@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Balanced',true,@(x) islogical(x) && isscalar(x));
p.addParameter('UseWeights',true,@(x) islogical(x) && isscalar(x));
p.addParameter('WeightComponent',0.60,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('LocalDistanceComponent',0.25,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('GlobalDistanceComponent',0.15,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('ADGlobalFloor',0.25,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('ReliabilityFloor',0.05,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.parse(varargin{:});
opt = p.Results;
K = opt.K;

neighborsIDs   = NaN(1,K);
neighborsSrc   = cell(1,K);
neighborsScore = NaN(1,K);
neighborsDC    = NaN(1,K);
neighborsW     = NaN(1,K);

if ~isfield(res,'CATMoS_winning_model_names') || numel(res.CATMoS_winning_model_names) < i || isempty(res.CATMoS_winning_model_names{i,1})
    error('res.CATMoS_winning_model_names{i,1} is missing or empty. Add/store winning model names in woe_corr first.');
end

winModels = res.CATMoS_winning_model_names{i,1};
if ischar(winModels) || isstring(winModels)
    winModels = cellstr(winModels);
end

modelNames = {'VT','NT','EPA','GHS','LD50'};
predModels = {predVT,predNT,predEPA,predGHS,predLD50};

AD_global = [get_res_value(res,'AD_VT',i), ...
             get_res_value(res,'AD_NT',i), ...
             get_res_value(res,'AD_EPA',i), ...
             get_res_value(res,'AD_GHS',i), ...
             get_res_value(res,'AD_LD50',i)];

AD_index = [get_res_value(res,'AD_index_VT',i), ...
            get_res_value(res,'AD_index_NT',i), ...
            get_res_value(res,'AD_index_EPA',i), ...
            get_res_value(res,'AD_index_GHS',i), ...
            get_res_value(res,'AD_index_LD50',i)];

Conf_index = [get_res_value(res,'Conf_index_VT',i), ...
              get_res_value(res,'Conf_index_NT',i), ...
              get_res_value(res,'Conf_index_EPA',i), ...
              get_res_value(res,'Conf_index_GHS',i), ...
              get_res_value(res,'Conf_index_LD50',i)];

AD_index   = clip01(AD_index);
Conf_index = clip01(Conf_index);

softAD = opt.ADGlobalFloor + (1-opt.ADGlobalFloor).*(AD_global == 1);
modelReliability = softAD .* AD_index .* Conf_index;
modelReliability(~isfinite(modelReliability)) = 0;
modelReliability = min(max(modelReliability,opt.ReliabilityFloor),1);

% Count maximum possible candidates from all winning endpoint models.
maxCandidates = 0;
for m = 1:numel(modelNames)
    if ismember(modelNames{m},winModels)
        ids = get_row(predModels{m}.neighbors,i);
        maxCandidates = maxCandidates + numel(ids);
    end
end

if maxCandidates == 0
    return;
end

candidateIDs    = NaN(1,maxCandidates);
candidateDC     = NaN(1,maxCandidates);
candidateW      = NaN(1,maxCandidates);
candidateScore  = NaN(1,maxCandidates);
candidateModel  = NaN(1,maxCandidates);
candidateSource = cell(1,maxCandidates);

pos = 0;

for m = 1:numel(modelNames)
    if ~ismember(modelNames{m},winModels)
        continue;
    end

    Pm = predModels{m};

    ids = get_row(Pm.neighbors,i);
    if isempty(ids)
        continue;
    end

    dc = get_neighbor_distances(Pm,i,ids);
    w  = get_neighbor_weights(Pm,i,numel(ids));

    localScore = local_neighbor_score(Pm,i,dc,w,opt);
    score = localScore .* modelReliability(m);

    n = numel(ids);
    idx = pos + (1:n);

    candidateIDs(idx)    = ids;
    candidateDC(idx)     = dc;
    candidateW(idx)      = w;
    candidateScore(idx)  = score;
    candidateModel(idx)  = m;
    candidateSource(idx) = repmat(modelNames(m),1,n);

    pos = pos + n;
end

candidateIDs    = candidateIDs(1:pos);
candidateDC     = candidateDC(1:pos);
candidateW      = candidateW(1:pos);
candidateScore  = candidateScore(1:pos);
candidateModel  = candidateModel(1:pos);
candidateSource = candidateSource(1:pos);

valid = isfinite(candidateIDs) & isfinite(candidateDC) & isfinite(candidateScore);

candidateIDs    = candidateIDs(valid);
candidateDC     = candidateDC(valid);
candidateW      = candidateW(valid);
candidateScore  = candidateScore(valid);
candidateModel  = candidateModel(valid);
candidateSource = candidateSource(valid);

if isempty(candidateIDs)
    return;
end

selectedCandIdx = NaN(1,K);
nSelected = 0;

% -------------------------------------------------------------------------
% Step 1: balanced first pass.
% Select one best unique neighbor from each winning endpoint model, ordered by
% endpoint reliability. This prevents all 5 neighbors from coming from one model.
% -------------------------------------------------------------------------
if opt.Balanced
    winningModelIdx = unique(candidateModel,'stable');
    relForWinning = modelReliability(winningModelIdx);
    [~,ordModels] = sort(relForWinning,'descend');
    winningModelIdx = winningModelIdx(ordModels);

    for jj = 1:numel(winningModelIdx)
        m = winningModelIdx(jj);

        cand = find(candidateModel == m);
        [~,ord] = sort(candidateScore(cand),'descend');
        cand = cand(ord);

        for kk = 1:numel(cand)
            id = candidateIDs(cand(kk));

            if ~ismember(id,candidateIDs(selectedCandIdx(1:nSelected)))
                nSelected = nSelected + 1;
                selectedCandIdx(nSelected) = cand(kk);
                break;
            end
        end

        if nSelected >= K
            break;
        end
    end
end

% -------------------------------------------------------------------------
% Step 2: fill remaining slots by final score.
% -------------------------------------------------------------------------
[~,ordAll] = sort(candidateScore,'descend');

for jj = 1:numel(ordAll)
    if nSelected >= K
        break;
    end

    c = ordAll(jj);
    id = candidateIDs(c);

    if ~ismember(id,candidateIDs(selectedCandIdx(1:nSelected)))
        nSelected = nSelected + 1;
        selectedCandIdx(nSelected) = c;
    end
end

% -------------------------------------------------------------------------
% Build outputs. If the same selected ID appeared in multiple models, source
% reports all sources, but score/DC/W come from the selected occurrence.
% -------------------------------------------------------------------------
for s = 1:nSelected
    c = selectedCandIdx(s);
    id = candidateIDs(c);
    sameID = candidateIDs == id;

    neighborsIDs(s)   = id;
    neighborsScore(s) = candidateScore(c);
    neighborsDC(s)    = candidateDC(c);
    neighborsW(s)     = candidateW(c);

    src = unique(candidateSource(sameID),'stable');
    neighborsSrc{s} = strjoin(src,'+');
end

end

% =========================================================================
% Helper functions
% =========================================================================

function row = get_row(A,i)
if isempty(A)
    row = [];
elseif isvector(A)
    row = A(:)';
else
    row = A(i,:);
end
end

function v = get_res_value(res,fieldName,i)
if isfield(res,fieldName) && numel(res.(fieldName)) >= i
    tmp = res.(fieldName);
    v = tmp(i,1);
else
    v = NaN;
end
end

function y = clip01(x)
y = x;
y(~isfinite(y)) = 0;
y = min(max(y,0),1);
end

function dc = get_neighbor_distances(Pm,i,ids)
% Prefer predX.dc if available. Otherwise extract from predX.D using ids.
if isfield(Pm,'dc') && ~isempty(Pm.dc)
    dc = get_row(Pm.dc,i);
else
    if ~isfield(Pm,'D') || isempty(Pm.D)
        dc = NaN(size(ids));
        return;
    end
    Drow = get_row(Pm.D,i);
    ids_int = round(ids);
    valid = ids_int >= 1 & ids_int <= numel(Drow) & isfinite(ids_int);
    dc = NaN(size(ids));
    dc(valid) = Drow(ids_int(valid));
end

if numel(dc) ~= numel(ids)
    dc = dc(:)';
    if numel(dc) > numel(ids)
        dc = dc(1:numel(ids));
    else
        dc = [dc NaN(1,numel(ids)-numel(dc))]; %#ok<AGROW>
    end
end
end

function w = get_neighbor_weights(Pm,i,n)
% Return raw source-model nearest-neighbor weights. These are used as the new
% neighborsW output and as one component of localScore.
if isfield(Pm,'w') && ~isempty(Pm.w)
    w = get_row(Pm.w,i);
else
    w = NaN(1,n);
end

w = w(:)';
if numel(w) ~= n
    if numel(w) > n
        w = w(1:n);
    else
        w = [w NaN(1,n-numel(w))]; %#ok<AGROW>
    end
end
end

function localScore = local_neighbor_score(Pm,i,dc,w,opt)
% Build a within-endpoint neighbor score in [0,1].
% Components:
%   1. normalized predX.w among the 5 neighbors
%   2. normalized local top-neighbor distance score
%   3. global distance percentile score from predX.D row, if available

n = numel(dc);
localScore = zeros(1,n);

% Normalize component weights.
comp = [opt.WeightComponent, opt.LocalDistanceComponent, opt.GlobalDistanceComponent];
if sum(comp) <= 0
    comp = [1,0,0];
end
comp = comp ./ sum(comp);

% --- Component 1: source-model weights ---
weightScore = zeros(1,n);
if opt.UseWeights && ~all(~isfinite(w))
    ww = w;
    ww(~isfinite(ww)) = 0;
    ww = max(ww,0);
    if max(ww) > 0
        weightScore = ww ./ max(ww);
    end
end

% If weights are unavailable/uninformative, give equal weight score rather than zero.
if max(weightScore) == 0
    weightScore = ones(1,n);
end

% --- Component 2: local distance rank/scale among selected neighbors ---
localDistScore = zeros(1,n);
validDC = isfinite(dc);
if any(validDC)
    dcv = dc;
    maxFinite = max(dcv(validDC));
    dcv(~validDC) = maxFinite;

    dmin = min(dcv);
    dmax = max(dcv);
    if dmax > dmin
        localDistScore = 1 - (dcv - dmin) ./ (dmax - dmin);
    else
        localDistScore = ones(1,n);
    end
end

% --- Component 3: global distance percentile from full distance matrix ---
globalDistScore = zeros(1,n);
if isfield(Pm,'D') && ~isempty(Pm.D)
    Drow = get_row(Pm.D,i);
    Drow = Drow(isfinite(Drow));
    if ~isempty(Drow)
        for k = 1:n
            if isfinite(dc(k))
                % Fraction of training compounds farther than or equal to this distance.
                % Higher score = closer relative to the endpoint's full distance distribution.
                globalDistScore(k) = mean(Drow >= dc(k));
            end
        end
    end
end

% If global distance component is unavailable, fold its mass into local distance.
if max(globalDistScore) == 0
    comp(2) = comp(2) + comp(3);
    comp(3) = 0;
    comp = comp ./ sum(comp);
end

localScore = comp(1).*weightScore + comp(2).*localDistScore + comp(3).*globalDistScore;
localScore(~isfinite(localScore)) = 0;
localScore = min(max(localScore,0),1);
end
