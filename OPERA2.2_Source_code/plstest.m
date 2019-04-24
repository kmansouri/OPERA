function pred = plstest(x,model)

% apply pls model on unknow samples (based on Frans Van Den Berg mypls routine)
%
% test = plstest(x,model)
%
% x         data matrix (n x p)
% model     pls model (see plsfit.m)
%
% test      predicted response (n x g)
%
% version 1.0 - september 2009
% Davide Ballabio
% Milano Chemometrics and QSAR Research Group
% www.disat.unimib.it/chm

W = model.W;
Q = model.Q;
P = model.P;
nF = size(model.T,2);
xscal = test_pretreatment(x,model.set.px);

yscal_c = 0;
for k = 1:nF
    Tnew = xscal*W(:,k)/(W(:,k)'*W(:,k));
    yscal_c = yscal_c + Tnew*Q(:,k)';
    xscal = xscal - Tnew*P(:,k)';
end
pred.yc = redo_scaling(yscal_c,model.set.py);

% % leverages
% xscal = test_pretreatment(x,model.set.px);
% T = model.T;
% Ttest = xscal*P(:,1:nF);
% pred.H = diag(Ttest*pinv(T'*T)*Ttest');
% pred.T = Ttest;