function [X_scal,param] = data_pretreatment(X,pret_type)

% pretreatment for multivariate data
%
% INPUT
% X:            data matrix (samples x variables)
% pret_type:    'cent' cenering
%               'scal' variance scaling
%               'auto' for autoscaling (centering + variance scaling)
%               'rang' range scaling (0-1)
%               'none' no scaling
%
% OUTPUT
% X_scal:       pretreated data matrix (samples x variables)
%
% version 1.0 - september 2009
% Davide Ballabio
% Milano Chemometrics and QSAR Research Group
% www.disat.unimib.it/chm


a = mean(X);
s = std(X);
m = min(X);
M = max(X);

if strcmp(pret_type,'cent')
    amat = repmat(a,size(X,1),1);
    X_scal = X - amat;
elseif strcmp(pret_type,'scal')
    f = find(s>0);
    smat = repmat(s,size(X,1),1);
    X_scal = zeros(size(X,1),size(X,2));
    X_scal = X(:,f)./smat(:,f);
elseif strcmp(pret_type,'auto')
    f = find(s>0);
    amat = repmat(a,size(X,1),1);
    smat = repmat(s,size(X,1),1);
    X_scal = zeros(size(X,1),size(X,2));
    X_scal(:,f) = (X(:,f) - amat(:,f))./smat(:,f);
elseif strcmp(pret_type,'rang')
    f = find(M - m > 0);
    mmat = repmat(m,size(X,1),1);
    Mmat = repmat(M,size(X,1),1);
    X_scal = zeros(size(X,1),size(X,2));
    X_scal(:,f) = (X(:,f) - mmat(:,f))./(Mmat(:,f) - mmat(:,f));       
else
    X_scal = X;
end

param.a = a;
param.s = s;
param.m = m;
param.M = M;
param.pret_type = pret_type;