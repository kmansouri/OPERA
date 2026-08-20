function res = woe_corrV76(res,i)
% =========================================================================
% woe_corr V7 : V5 cascade with a probabilistic tie resolver at the bottom.
%
% Identical to V5 (woe_corr_full) for:
%   §1-§6  matrix construction (votes, weights, AD, ADi, Cfi)
%   §7-overrides (extreme, cat2, nontoxic, bin6)
%   §7-AD  all-5-in-domain gate
%   §7a    LD50 gated upward relax
%   §7b    low-bin safety rule
%   §7c    high-bin safety rule
%   §8     LD50 re-centering with clamp into final bin
%   §9     combined AD / confidence outputs
%
% ONLY DIFFERENCE from V5: the §7d fallback tie resolver.
%   V5 §7d: per-model reliability via AD*(AD_index+Conf)/2, summed per bin.
%   V7 §7d: per-model SOFT DISTRIBUTION restricted to the tied bins, weighted
%           by MULTIPLICATIVE reliability (AD_global*AD_index*Conf), pooled
%           linearly. Pick the argmax of the pooled posterior over the tied
%           set. Falls back to min(maxbins) if the pooled posterior is empty.
%
% Why this is a strictly-better fallback than V5's:
%   * Soft per-model distributions respect ordinal proximity: LD50's point
%     estimate gives more mass to the closer tied bin than the farther one,
%     even if both are inside its +/-0.25 band. V5's binary vote treated them
%     as equal.
%   * Multiplicative reliability (AD_global * AD_index * Conf) requires BOTH
%     applicability AND fidelity; additive let one mask the other.
%   * No new heuristics added — the existing safety/LD50 rules (§7a-§7c) run
%     first, so the probabilistic resolver only handles what's left over.
% =========================================================================

%% ------------------ V7 §7d tunables --------------------------------------
v7_LD50_sigma     = 0.25;   % LD50 kernel std (= design uncertainty band)
v7_cat_tau_base   = 0.45;   % base spread of categorical kernels (bins)
v7_cat_tau_unconf = 0.90;   % extra spread = (1-conf)*this for categoricals
v7_bin_tau        = 1.10;   % spread of binary VT/NT kernels
v7_adGlobalFloor  = 0.25;   % out-of-domain models keep this fraction of weight

%% ------------------ Optional behaviour toggles (V5 defaults) -------------
USE_GATED_LD50_TIEBREAK = true;
USE_WEIGHTED_NEAR_TIE   = true;
USE_LD50_BIN_CLAMP      = true;
USE_ABSTENTION_FLAG     = false;
ABSTAIN_AD_THRESH       = 0.30;

% Last-edit safety controls.
% - Cat2 override: VT alone no longer forces WoE=2 unless bin 2 has
%   stronger vote support.
% - LD50 hard overrides/tie-breaks require LD50 to be reliable.
% - Endpoint Conf_index values are bounded internally, so final
%   Conf_index_CATMoS cannot exceed 1 because of upstream Conf_index >1.
USE_CONF_CLAMP            = true;
USE_PROB_CONF_TERM        = true;
LD50_RELIABLE_ADI_MIN     = 0.50;
LD50_RELIABLE_CONF_MIN    = 0.50;

% LD50 downward-anchor controls. These are intentionally weaker than the
% full LD50 reliability gate and are used only to prevent upward/non-toxic
% drift when LD50's own bin/span is lower than competing near-tie bins.
% They do NOT allow a weak high LD50 to push the consensus upward.
USE_LD50_DOWNWARD_ANCHOR = true;
LD50_ANCHOR_ADI_MIN      = 0.20;
LD50_ANCHOR_CONF_MIN     = 0.25;

% Extra downward-only protection against upward/non-toxic drift. These rules
% only act when the final/tied candidate is higher than the LD50 +/-0.25 span
% and the LD50-supported anchor bin has enough vote/category support. They
% do NOT let a weak high LD50 push the consensus upward.
USE_LD50_SPAN_OVERLAP_ANCHOR = true;
USE_LD50_UPWARD_JUMP_GUARD   = true;
LD50_UPWARD_GUARD_MAX_VOTE_GAP = 1;    % anchor bin must be within this many votes of high bin
LD50_UPWARD_GUARD_MIN_CAT_SUPPORT = 1;  % categorical rows supporting anchor bin

%% =========================== §1  INITIALISE ==============================
M   = zeros(5,7); W = zeros(5,7); AD = zeros(5,7); ADi = zeros(5,7); Cfi = zeros(5,7);

% Bounded confidence copies used only inside this consensus function. The raw
% res.Conf_index_* fields are not overwritten.
if USE_CONF_CLAMP
    confVT   = clip01(res.Conf_index_VT(i,1));
    confNT   = clip01(res.Conf_index_NT(i,1));
    confEPA  = clip01(res.Conf_index_EPA(i,1));
    confGHS  = clip01(res.Conf_index_GHS(i,1));
    confLD50 = clip01(res.Conf_index_LD50(i,1));
else
    confVT   = res.Conf_index_VT(i,1);
    confNT   = res.Conf_index_NT(i,1);
    confEPA  = res.Conf_index_EPA(i,1);
    confGHS  = res.Conf_index_GHS(i,1);
    confLD50 = res.Conf_index_LD50(i,1);
end

%% =========================== §2  VT votes ================================
if res.CATMoS_VT_pred(i,1)==1
    M(1,1:2)=ones(1,2);
    W(1,1:2)=repelem((res.AD_index_VT(i,1)+confVT)*res.AD_VT(i,1)/2,2);
    AD(1,1:2)=repelem(res.AD_VT(i,1),2); ADi(1,1:2)=repelem(res.AD_index_VT(i,1),2); Cfi(1,1:2)=repelem(confVT,2);
else
    M(1,3:end)=ones(1,5);
    W(1,3:end)=repelem((res.AD_index_VT(i,1)+confVT)*res.AD_VT(i,1)/2,5);
    AD(1,3:end)=repelem(res.AD_VT(i,1),5); ADi(1,3:end)=repelem(res.AD_index_VT(i,1),5); Cfi(1,3:end)=repelem(confVT,5);
end

%% =========================== §3  NT votes ================================
if res.CATMoS_NT_pred(i,1)==1
    M(2,6:7)=ones(1,2);
    W(2,6:7)=repelem((res.AD_index_NT(i,1)+confNT)*res.AD_NT(i,1)/2,2);
    AD(2,6:7)=repelem(res.AD_NT(i,1),2); ADi(2,6:7)=repelem(res.AD_index_NT(i,1),2); Cfi(2,6:7)=repelem(confNT,2);
else
    M(2,1:5)=ones(1,5);
    W(2,1:5)=repelem((res.AD_index_NT(i,1)+confNT)*res.AD_NT(i,1)/2,5);
    AD(2,1:5)=repelem(res.AD_NT(i,1),5); ADi(2,1:5)=repelem(res.AD_index_NT(i,1),5); Cfi(2,1:5)=repelem(confNT,5);
end

%% =========================== §4  EPA votes ===============================
switch(res.CATMoS_EPA_pred(i,1))
    case 1, cols=1:2; case 2, cols=3:4; case 3, cols=5:6; case 4, cols=7; otherwise, cols=[];
end
if ~isempty(cols)
    n=numel(cols); M(3,cols)=ones(1,n);
    W(3,cols)=repelem((res.AD_index_EPA(i,1)+confEPA)*res.AD_EPA(i,1)/2,n);
    AD(3,cols)=repelem(res.AD_EPA(i,1),n); ADi(3,cols)=repelem(res.AD_index_EPA(i,1),n); Cfi(3,cols)=repelem(confEPA,n);
end

%% =========================== §5  GHS votes ===============================
switch(res.CATMoS_GHS_pred(i,1))
    case 1, cols=1; case 2, cols=2; case 3, cols=3; case 4, cols=4:5; case 5, cols=6:7; otherwise, cols=[];
end
if ~isempty(cols)
    n=numel(cols); M(4,cols)=ones(1,n);
    W(4,cols)=repelem((res.AD_index_GHS(i,1)+confGHS)*res.AD_GHS(i,1)/2,n);
    AD(4,cols)=repelem(res.AD_GHS(i,1),n); ADi(4,cols)=repelem(res.AD_index_GHS(i,1),n); Cfi(4,cols)=repelem(confGHS,n);
end

%% =========================== §6  LD50 votes ==============================
if ~isnan(res.CATMoS_LD50_pred(i,1))
    LD50_edges = log10([5,50,300,500,2000,5000]);
    A = find((res.CATMoS_LD50_pred(i,1)-0.25) < LD50_edges, 1, 'first');
    B = find((res.CATMoS_LD50_pred(i,1)+0.25) < LD50_edges, 1, 'first');
    if isempty(A), A = 7; end
    if isempty(B), B = 7; end
    n=B-A+1;
    M(5,A:B)=ones(1,n);
    W(5,A:B)=repelem((res.AD_index_LD50(i,1)+confLD50)*res.AD_LD50(i,1)/2,n);
    AD(5,A:B)=repelem(res.AD_LD50(i,1),n); ADi(5,A:B)=repelem(res.AD_index_LD50(i,1),n); Cfi(5,A:B)=repelem(confLD50,n);
end

%% ===================== §7  CONSENSUS DECISION CASCADE ====================
votes   = sum(M);
maxbins = find(votes==max(votes));
colGe2  = find(votes>=2);

% LD50 hard upward/non-toxic rules are allowed only when LD50 is reliable.
% LD50 still votes in M and contributes softly in §7d/§8b when this is false.
LD50_reliable = ~isnan(res.CATMoS_LD50_pred(i,1)) && ...
    res.AD_LD50(i,1)==1 && ...
    res.AD_index_LD50(i,1)>=LD50_RELIABLE_ADI_MIN && ...
    confLD50>=LD50_RELIABLE_CONF_MIN;

% Central LD50 bin and LD50 +/-0.25 span from the LD50 vote row.
LD50_bin = [];
LD50_span = [];
if ~isnan(res.CATMoS_LD50_pred(i,1))
    LD50_bin = find((res.CATMoS_LD50_pred(i,1)-log10([5,50,300,500,2000,5000,Inf]))<0, 1);
    LD50_span = find(M(5,:) ~= 0);
end

% Less strict LD50 use: allowed only as a downward/conservative anchor.
% This prevents a bin-5-ish LD50 from being pushed to bin 6/7 by NT/high-bin
% categorical votes, while still blocking an unreliable high LD50 from
% pushing the consensus upward/non-toxic.
LD50_anchor_available = USE_LD50_DOWNWARD_ANCHOR && ...
    ~isempty(LD50_bin) && ...
    res.AD_index_LD50(i,1)>=LD50_ANCHOR_ADI_MIN && ...
    confLD50>=LD50_ANCHOR_CONF_MIN;

% Cat2 override components. Direct severe evidence can still trigger bin 2,
% but VT alone now needs at least 3 bin-2 votes, preventing VT+NT alone from
% forcing an over-toxic consensus in broad conflicts.
cat2_direct_evidence = ...
       res.CATMoS_GHS_pred(i,1)<=2 ...
    || (LD50_reliable && res.CATMoS_LD50_pred(i,1)<=log10(75)) ...
    || res.CATMoS_EPA_pred(i,1)==1;

cat2_vt_supported = ...
       res.CATMoS_VT_pred(i,1)==1 ...
    && confVT>=0.55 ...
    && votes(2)>=3;

% ---- §7-overrides ----
if   ( votes(1)>=3 && min(maxbins)<=3 && ((res.CATMoS_GHS_pred(i,1)==1 && confGHS>=0.5) || (LD50_reliable && res.CATMoS_LD50_pred(i,1)<=log10(17.5))) ) ...
  || ( votes(2)>=3 && min(maxbins)<=3 && ((res.CATMoS_GHS_pred(i,1)==1 && confGHS>=0.5) || (LD50_reliable && res.CATMoS_LD50_pred(i,1)<=log10(17.5))) )
    WoE=1;
elseif ( votes(2)>=4 && min(maxbins)<=3 ) ...
    || ( votes(2)>=2 && min(maxbins)<=3 && (cat2_direct_evidence || cat2_vt_supported) )
    WoE=2;
elseif votes(7)>=3 && LD50_reliable && res.CATMoS_LD50_pred(i,1)>=log10(4000) && confLD50>=confEPA
    WoE=7;
elseif ( votes(6)>=4 && max(maxbins)>=6 ) && res.CATMoS_NT_pred(i,1)==1 && LD50_reliable && res.CATMoS_LD50_pred(i,1)>=log10(2500) && res.CATMoS_LD50_pred(i,1)<=log10(5000) && res.CATMoS_EPA_pred(i,1)==3
    WoE=6;

% ---- §7-AD: all-5-in-domain gate ----
elseif (res.AD_VT(i,1)==1 && res.AD_NT(i,1)==1 && res.AD_EPA(i,1)==1 && ...
        res.AD_GHS(i,1)==1 && res.AD_LD50(i,1)==1) && ...
        res.CATMoS_LD50_pred(i,1)<=log10(500) && ~isempty(colGe2) && min(colGe2)<=4
    WoE=find([ res.CATMoS_LD50_pred(i,1)<=log10(5) ...
               res.CATMoS_LD50_pred(i,1)>log10(5)    && res.CATMoS_LD50_pred(i,1)<=log10(50) ...
               res.CATMoS_LD50_pred(i,1)>log10(50)   && res.CATMoS_LD50_pred(i,1)<=log10(300) ...
               res.CATMoS_LD50_pred(i,1)>log10(300)  && res.CATMoS_LD50_pred(i,1)<=log10(500) ...
               res.CATMoS_LD50_pred(i,1)>log10(500)  && res.CATMoS_LD50_pred(i,1)<=log10(2000) ...
               res.CATMoS_LD50_pred(i,1)>log10(2000) && res.CATMoS_LD50_pred(i,1)<=log10(5000) ]);
    if WoE > min(colGe2), WoE=min(colGe2); end

% ---- §7-plurality + tie-breaking ----
else
    if USE_WEIGHTED_NEAR_TIE
        WoE = find(votes >= max(votes)-1);   % near-ties (within 1 vote)
    else
        WoE = maxbins;
    end

    if length(WoE)>1
        % §7a  LD50 tie-break.
        % Reliable LD50 can tie-break normally. Less-reliable LD50 can only
        % act as a downward anchor when its central bin is one of the
        % candidate bins and is lower than the upper tied/near-tied bins.
        LD50_tiebreak_allowed = ...
               ~isempty(LD50_bin) ...
            && ismember(LD50_bin,WoE) ...
            && (LD50_reliable || (LD50_anchor_available && LD50_bin < max(WoE)));

        % §7c high-bin blocker, computed before the if/elseif chain below.
        % It blocks only upward drift to bin 6/7 when LD50's own +/-0.25
        % span is lower and overlaps the tied/near-tied candidates.
        LD50_blocks_high_tie = ...
               LD50_anchor_available ...
            && ~isempty(LD50_span) ...
            && any(ismember(WoE,LD50_span)) ...
            && max(LD50_span) < 6 ...
            && max(WoE) >= 6;

        % §7a2 candidate: if LD50 +/-0.25 overlaps a lower tied/near-tied
        % bin, keep that highest LD50-supported overlap as a downward anchor.
        % This catches cases that the central-bin tiebreak misses and avoids
        % falling through to §7d where categorical high-bin evidence may win.
        LD50_span_overlap_anchor_bin = [];
        if USE_LD50_SPAN_OVERLAP_ANCHOR && LD50_anchor_available && ...
                ~isempty(LD50_span) && max(WoE) > max(LD50_span)
            overlapBins = intersect(WoE,LD50_span);
            if ~isempty(overlapBins)
                candidateAnchor = max(overlapBins);
                if candidateAnchor < max(WoE) && ...
                        sum(M(1:4,candidateAnchor)) >= LD50_UPWARD_GUARD_MIN_CAT_SUPPORT && ...
                        votes(candidateAnchor) >= max(votes(WoE)) - LD50_UPWARD_GUARD_MAX_VOTE_GAP
                    LD50_span_overlap_anchor_bin = candidateAnchor;
                end
            end
        end

        if LD50_tiebreak_allowed
            if ~LD50_reliable && LD50_bin < max(WoE)
                WoE = LD50_bin;
            elseif USE_GATED_LD50_TIEBREAK
                relax = [];
                for b = WoE(:)'
                    if b>LD50_bin && ismember(b,LD50_span)
                        catRows = find(M(1:4,b));
                        if numel(catRows)>=2 && mean(Cfi(catRows,b))>=confLD50
                            relax(end+1)=b; %#ok<AGROW>
                        end
                    end
                end
                if ~isempty(relax)
                    v = arrayfun(@(b) numel(find(M(1:4,b))), relax);
                    relax = relax(v==max(v));
                    WoE = min(relax);
                else
                    WoE = LD50_bin;
                end
            else
                WoE = LD50_bin;
            end

        % §7a2 LD50-span overlap downward anchor.
        elseif ~isempty(LD50_span_overlap_anchor_bin)
            WoE = LD50_span_overlap_anchor_bin;

        % §7b  low-bin tie
        elseif min(WoE)<=3 && max(WoE)<=4 && (LD50_reliable || LD50_anchor_available) && res.CATMoS_LD50_pred(i,1)<=log10(500)
            if min(WoE)>1
                WoE=min(WoE);
            elseif res.CATMoS_GHS_pred(i,1)==1 || res.CATMoS_LD50_pred(i,1)<=log10(100)
                WoE=min(WoE);
            elseif res.CATMoS_LD50_pred(i,1)<=log10(300)
                WoE=round(mean(WoE));
            else
                WoE=max(WoE);
            end

        % §7c  high-bin tie
        % Do not let the high-bin rule jump to bin 6/7 when LD50's own
        % +/-0.25 span is lower, overlaps the candidate set, and does not
        % support bin 6. This fixes bin-5-ish LD50 cases being recentered
        % upward into bin 6 solely because NT/high-bin votes were present.
        elseif max(WoE)>=6 && ~LD50_blocks_high_tie && ...
               (res.CATMoS_NT_pred(i,1)==1 || (LD50_reliable && res.CATMoS_LD50_pred(i,1)>=log10(2500)))
            if ismember(6,WoE), WoE=6; else, WoE=max(WoE); end

        % §7d  PROBABILISTIC tie resolver  (V7 only -- replaces V5's W-weighted pick)
        %
        % Build a per-model soft distribution over the tied bins, weight by
        % multiplicative reliability, pool linearly, pick the argmax.
        else
            tied = WoE(:)';
            nTied = numel(tied);

            % per-model multiplicative reliability with soft floor for out-of-domain
            softg = @(g) v7_adGlobalFloor + (1-v7_adGlobalFloor)*(g==1);
            relW = [ softg(res.AD_VT(i,1))  * res.AD_index_VT(i,1)  * confVT  , ...
                     softg(res.AD_NT(i,1))  * res.AD_index_NT(i,1)  * confNT  , ...
                     softg(res.AD_EPA(i,1)) * res.AD_index_EPA(i,1) * confEPA , ...
                     softg(res.AD_GHS(i,1)) * res.AD_index_GHS(i,1) * confGHS , ...
                     softg(res.AD_LD50(i,1))* res.AD_index_LD50(i,1)* confLD50 ];

            % per-model soft distribution evaluated on the tied bins only
            Dt = zeros(5,nTied);

            % VT / NT — coarse half-line kernels
            VTset = ternarySet(res.CATMoS_VT_pred(i,1)==1,[1 2],[3 4 5 6 7]);
            NTset = ternarySet(res.CATMoS_NT_pred(i,1)==1,[6 7],[1 2 3 4 5]);
            Dt(1,:) = local_kernel(tied,VTset,v7_bin_tau);
            Dt(2,:) = local_kernel(tied,NTset,v7_bin_tau);

            % EPA / GHS — categorical kernels, spread widens as conf drops
            switch res.CATMoS_EPA_pred(i,1)
                case 1, EPAset=[1 2]; case 2, EPAset=[3 4]; case 3, EPAset=[5 6]; case 4, EPAset=7; otherwise, EPAset=[];
            end
            switch res.CATMoS_GHS_pred(i,1)
                case 1, GHSset=1; case 2, GHSset=2; case 3, GHSset=3; case 4, GHSset=[4 5]; case 5, GHSset=[6 7]; otherwise, GHSset=[];
            end
            Dt(3,:) = local_kernel(tied,EPAset,v7_cat_tau_base + v7_cat_tau_unconf*(1-confEPA));
            Dt(4,:) = local_kernel(tied,GHSset,v7_cat_tau_base + v7_cat_tau_unconf*(1-confGHS));

            % LD50 — true Gaussian mass per bin from its calibrated +/-sigma band
            if ~isnan(res.CATMoS_LD50_pred(i,1))
                edgesLo = [-inf, log10([5,50,300,500,2000,5000])];
                edgesHi = [log10([5,50,300,500,2000,5000]), inf];
                full = 0.5*(1+erf((edgesHi - res.CATMoS_LD50_pred(i,1))/(v7_LD50_sigma*sqrt(2)))) ...
                     - 0.5*(1+erf((edgesLo - res.CATMoS_LD50_pred(i,1))/(v7_LD50_sigma*sqrt(2))));
                if sum(full)>0, full = full/sum(full); else, full = ones(1,7)/7; end
                Dt(5,:) = full(tied);
                if sum(Dt(5,:))>0, Dt(5,:) = Dt(5,:)/sum(Dt(5,:)); end
            else
                Dt(5,:) = 0; relW(5) = 0;
            end

            % pool linearly with reliability weights
            if sum(relW)<=0, relW = ones(1,5); end
            pooled = relW * Dt;                          % 1 x nTied

            if sum(pooled)>0
                [~,kbest] = max(pooled);
                WoE = tied(kbest);
            else
                WoE = min(maxbins);                       % defensive fallback
            end
        end

    else
        if WoE==4 && votes(4)-votes(3)<=1 && (LD50_reliable || LD50_anchor_available) && res.CATMoS_LD50_pred(i,1)<=log10(500)
            WoE=WoE-1;
        end
    end
end

% Final downward-only guard for scalar outcomes. This catches cases where
% the cascade/fallback still selected bin 6 even though LD50 +/-0.25 only
% supports lower bins and the LD50-supported anchor bin is close in votes.
% It is deliberately narrow to avoid undoing strong high-bin consensus.
if USE_LD50_UPWARD_JUMP_GUARD && numel(WoE)==1 && ...
        ~isempty(LD50_span) && WoE >= 6 && max(LD50_span) < WoE
    anchorBin = max(LD50_span);
    if anchorBin >= 1 && anchorBin <= 7 && ...
            sum(M(1:4,anchorBin)) >= LD50_UPWARD_GUARD_MIN_CAT_SUPPORT && ...
            votes(anchorBin) >= votes(WoE) - LD50_UPWARD_GUARD_MAX_VOTE_GAP
        WoE = anchorBin;
    end
end

WoE = min(WoE);
%res.CATMoS_WoE(i,1) = WoE;

%% =============== §8b PROBABILISTIC POSTERIOR (for Conf_index) ============
% Build the full 7-bin posterior using the same soft-distribution machinery
% the §7d resolver uses, but evaluated over ALL bins. This is computed HERE
% (BEFORE §8 overwrites the categorical predictions) so it reflects the
% ORIGINAL model inputs, not the harmonized outputs. Used ONLY by
% Conf_index_CATMoS in §9 — does NOT affect the WoE choice in any way.

allBins = 1:7;
softg_f = @(g) v7_adGlobalFloor + (1-v7_adGlobalFloor)*(g==1);
relW_f = [ softg_f(res.AD_VT(i,1))  * res.AD_index_VT(i,1)  * confVT  , ...
           softg_f(res.AD_NT(i,1))  * res.AD_index_NT(i,1)  * confNT  , ...
           softg_f(res.AD_EPA(i,1)) * res.AD_index_EPA(i,1) * confEPA , ...
           softg_f(res.AD_GHS(i,1)) * res.AD_index_GHS(i,1) * confGHS , ...
           softg_f(res.AD_LD50(i,1))* res.AD_index_LD50(i,1)* confLD50 ];

Df = zeros(5,7);
Df(1,:) = local_kernel(allBins, ternarySet(res.CATMoS_VT_pred(i,1)==1,[1 2],[3 4 5 6 7]), v7_bin_tau);
Df(2,:) = local_kernel(allBins, ternarySet(res.CATMoS_NT_pred(i,1)==1,[6 7],[1 2 3 4 5]), v7_bin_tau);
switch res.CATMoS_EPA_pred(i,1)
    case 1, EPAset_f=[1 2]; case 2, EPAset_f=[3 4]; case 3, EPAset_f=[5 6]; case 4, EPAset_f=7; otherwise, EPAset_f=[];
end
switch res.CATMoS_GHS_pred(i,1)
    case 1, GHSset_f=1; case 2, GHSset_f=2; case 3, GHSset_f=3; case 4, GHSset_f=[4 5]; case 5, GHSset_f=[6 7]; otherwise, GHSset_f=[];
end
Df(3,:) = local_kernel(allBins, EPAset_f, v7_cat_tau_base + v7_cat_tau_unconf*(1-confEPA));
Df(4,:) = local_kernel(allBins, GHSset_f, v7_cat_tau_base + v7_cat_tau_unconf*(1-confGHS));

if ~isnan(res.CATMoS_LD50_pred(i,1))
    edgesLo_f = [-inf, log10([5,50,300,500,2000,5000])];
    edgesHi_f = [log10([5,50,300,500,2000,5000]), inf];
    LD_full = 0.5*(1+erf((edgesHi_f - res.CATMoS_LD50_pred(i,1))/(v7_LD50_sigma*sqrt(2)))) ...
            - 0.5*(1+erf((edgesLo_f - res.CATMoS_LD50_pred(i,1))/(v7_LD50_sigma*sqrt(2))));
    if sum(LD_full)>0, Df(5,:) = LD_full/sum(LD_full); else, Df(5,:) = ones(1,7)/7; end
else
    Df(5,:) = 0; relW_f(5) = 0;
end

if sum(relW_f)<=0, relW_f = ones(1,5); end
P_full = relW_f * Df;
if sum(P_full)>0, P_full = P_full / sum(P_full); else, P_full = ones(1,7)/7; end

%% =============== §8  HARMONISE PREDICTIONS + RE-CENTRE LD50 ==============
VT_by_WoE  = [1,1,0,0,0,0,0];
NT_by_WoE  = [0,0,0,0,0,1,1];
EPA_by_WoE = [1,1,2,2,3,3,4];
GHS_by_WoE = [1,2,3,4,4,5,5];
res.CATMoS_VT_pred(i,1) = VT_by_WoE(WoE);
res.CATMoS_NT_pred(i,1) = NT_by_WoE(WoE);
res.CATMoS_EPA_pred(i,1)= EPA_by_WoE(WoE);
res.CATMoS_GHS_pred(i,1)= GHS_by_WoE(WoE);

res.CATMoS_LD50_pred(i,1) = adjust_ld50_to_woe(res.CATMoS_LD50_pred(i,1), WoE);

if USE_LD50_BIN_CLAMP
    edges_lo = [ -Inf, log10([5 50 300 500 2000 5000]) ];
    edges_hi = [ log10([5 50 300 500 2000 5000]), Inf   ];
    res.CATMoS_LD50_pred(i,1) = min( max(res.CATMoS_LD50_pred(i,1), edges_lo(WoE)), edges_hi(WoE) );
end

%% =============== §9  COMBINED AD / CONFIDENCE (+ flag) ===================
rows = M(:,WoE) ~= 0;

% Winning endpoint models for the final consensus bin.
% Row mapping: 1=VT, 2=NT, 3=EPA, 4=GHS, 5=LD50.
modelNames = {'VT','NT','EPA','GHS','LD50'};
winningModelIdx = find(rows)';
%res.CATMoS_winning_model_mask(i,:) = false(1,5);
% if ~isempty(winningModelIdx)
%     res.CATMoS_winning_model_mask(i,winningModelIdx) = true;
% end
% res.CATMoS_winning_model_count(i,1) = numel(winningModelIdx);
res.CATMoS_winning_model_idx{i,1} = winningModelIdx;
res.CATMoS_winning_model_names{i,1} = modelNames(winningModelIdx);

if any(rows)
    res.AD_CATMoS(i,1)        = round(mean(AD(rows,WoE)));
    res.AD_index_CATMoS(i,1)  = max(mean(ADi(rows,WoE)),  median(ADi(rows,WoE)));
    % mass_within1: reliability-weighted soft mass within +/-1 bin of WoE.
    % Replaces the old hard-vote term sum(M(:,WoE))/5 as the agreement
    % component of Conf_index_CATMoS. Better because: (a) continuous not
    % discrete-5-step, (b) reliability-weighted, (c) ordinal-aware (adjacent
    % bins support the estimate), (d) flags fragile predictions the vote
    % count masked (mol2, mol14 correctly deflated; mol8 correctly kept high).
    mass_within1 = sum(P_full(max(WoE-1,1):min(WoE+1,7)));
    if isnan(mass_within1) || mass_within1 < 0 || mass_within1 > 1
        mass_within1 = sum(M(:,WoE))/5;   % defensive fallback to old term
    end
    confVals = Cfi(rows,WoE);
    if USE_CONF_CLAMP
        confVals = clip01(confVals);
    end
    hardAgreement = sum(M(:,WoE))/5;
    if USE_PROB_CONF_TERM
        agreementTerm = mass_within1;
    else
        agreementTerm = hardAgreement;
    end
    res.Conf_index_CATMoS(i,1)= median([mean(confVals), median(confVals), agreementTerm]);
    res.Conf_index_CATMoS(i,1)= min(max(res.Conf_index_CATMoS(i,1),0),1);
else
    res.AD_CATMoS(i,1)=NaN; res.AD_index_CATMoS(i,1)=NaN; res.Conf_index_CATMoS(i,1)=NaN;
end

if USE_ABSTENTION_FLAG
    meanADi = mean([res.AD_index_VT(i,1) res.AD_index_NT(i,1) res.AD_index_EPA(i,1) ...
                    res.AD_index_GHS(i,1) res.AD_index_LD50(i,1)]);
    res.CATMoS_unreliable(i,1) = double(meanADi < ABSTAIN_AD_THRESH);
end

end


%% ========================================================================
%  Helpers
%% ========================================================================
function s = ternarySet(cond,a,b)
if cond, s=a; else, s=b; end
end

function d = local_kernel(tiedBins,prefSet,tau)
% Soft distribution over tiedBins: mass decays with squared distance from
% the preferred bin-set. Empty set -> uniform.
if isempty(prefSet)
    d = ones(1,numel(tiedBins))/numel(tiedBins); return;
end
dist = arrayfun(@(b) min(abs(b-prefSet)), tiedBins);
d = exp(-(dist.^2)/(2*tau^2));
if sum(d)>0, d = d/sum(d); else, d = ones(1,numel(tiedBins))/numel(tiedBins); end
end

function ld50 = adjust_ld50_to_woe(ld50, WoE)
% Original CATMoS LD50 re-centering formulas (unchanged).
if isnan(ld50), return; end
l5=log10(5); l50=log10(50); l300=log10(300); l500=log10(500);
l2000=log10(2000); l5000=log10(5000); l10000=log10(10000);
switch WoE
    case 1
        if ld50 > l5
            if ld50-0.25 <= l5, ld50=(ld50-0.25+l5)/2;
            elseif (ld50-1)/0.65 <= l5, ld50=(ld50-1)/0.65;
            else, ld50=l5*ld50/4; end
        end
    case 2
        if ld50 > l50
            if ld50-0.25 <= l50 && (ld50-0.25+l50)/2 >= l5, ld50=(ld50-0.25+l50)/2;
            elseif (ld50-1)/0.65 <= l50 && (ld50-1)/0.65 >= l5, ld50=(ld50-1)/0.65;
            else, ld50=l5+(l50-l5)*ld50/4; end
        elseif ld50 < l5
            if ld50+0.25 >= l5 && (ld50+0.25+l5)/2 <= l50, ld50=(ld50+0.25+l5)/2;
            elseif (ld50-0.781)/0.7369 >= l5 && (ld50-0.781)/0.7369 <= l50, ld50=(ld50-0.781)/0.7369;
            else, ld50=l5+(l50-l5)*ld50/4; end
        end
    case 3
        if ld50 > l300
            if ld50-0.25 <= l300 && (ld50-0.25+l300)/2 >= l50, ld50=(ld50-0.25+l300)/2;
            elseif (ld50-0.781)/0.7369 <= l300 && (ld50-0.781)/0.7369 >= l50, ld50=(ld50-0.781)/0.7369;
            else, ld50=l50+(l300-l50)*ld50/4; end
        elseif ld50 < l50
            if ld50+0.25 >= l50 && (ld50+0.25+l50)/2 <= l300, ld50=(ld50+0.25+l50)/2;
            elseif (ld50-0.781)/0.7369 >= l50 && (ld50-0.781)/0.7369 <= l300, ld50=(ld50-0.781)/0.7369;
            else, ld50=l50+(l300-l50)*ld50/4; end
        end
    case 4
        if ld50 > l500
            if ld50-0.25 <= l500 && (ld50-0.25+l500)/2 >= l300, ld50=(ld50-0.25+l500)/2;
            elseif (ld50-0.781)/0.7369 <= l500 && (ld50-0.781)/0.7369 >= l300, ld50=(ld50-0.781)/0.7369;
            else, ld50=l300+(l500-l300)*ld50/4; end
        elseif ld50 < l300
            if ld50+0.25 >= l300 && (ld50+0.25+l300)/2 <= l500, ld50=(ld50+0.25+l300)/2;
            elseif (ld50-0.781)/0.7369 >= l300 && (ld50-0.781)/0.7369 <= l500, ld50=(ld50-0.781)/0.7369;
            else, ld50=l300+(l500-l300)*ld50/4; end
        end
    case 5
        if ld50 > l2000
            if ld50-0.25 <= l2000 && (ld50-0.25+l2000)/2 >= l500, ld50=(ld50-0.25+l2000)/2;
            elseif (ld50-0.781)/0.7369 <= l2000 && (ld50-0.781)/0.7369 >= l500, ld50=(ld50-0.781)/0.7369;
            else, ld50=l500+(l2000-l500)*ld50/4; end
        elseif ld50 < l500
            if ld50+0.25 >= l500 && (ld50+0.25+l500)/2 <= l2000, ld50=(ld50+0.25+l500)/2;
            elseif (ld50-0.781)/0.7369 >= l500 && (ld50-0.781)/0.7369 <= l2000, ld50=(ld50-0.781)/0.7369;
            else, ld50=l500+(l2000-l500)*ld50/4; end
        end
    case 6
        if ld50 > l5000 && (ld50-0.25+l5000)/2 >= l2000
            if ld50-0.25 <= l5000, ld50=(ld50-0.25+l5000)/2;
            elseif (ld50-0.781)/0.7369 <= l5000 && (ld50-0.781)/0.7369 >= l2000, ld50=(ld50-0.781)/0.7369;
            else, ld50=l2000+(l5000-l2000)*ld50/4; end
        elseif ld50 < l2000
            if ld50+0.25 >= l2000 && (ld50+0.25+l2000)/2 <= l5000, ld50=(ld50+0.25+l2000)/2;
            elseif (ld50-0.781)/0.7369 >= l2000 && (ld50-0.781)/0.7369 <= l5000, ld50=(ld50-0.781)/0.7369;
            else, ld50=l2000+(l5000-l2000)*ld50/4; end
        end
    case 7
        if ld50 < l5000
            if ld50+0.25 >= l5000, ld50=(ld50+0.25+l5000)/2;
            elseif (ld50-0.781)/0.7369 >= l5000, ld50=(ld50-0.781)/0.7369;
            else, ld50=l5000+(l10000-l5000)*ld50/4; end
        end
end
end

function y = clip01(x)
%CLIP01 Clamp numeric scalar/vector/matrix values to [0,1].
% Non-finite entries are treated as 0.
y = x;
y(~isfinite(y)) = 0;
y = min(max(y,0),1);
end
