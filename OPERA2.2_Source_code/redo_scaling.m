function X = redo_scaling(X_scal,param)

% redo scaling (from scaled data to original data)
%
% INPUT
% X_scal:   pretreated data matrix (samples x variables)
% param:    output data structure from data_pretreatment routine
%
% OUTPUT
% X:        data matrix (samples x variables)
%
% version 1.0 - september 2009
% Davide Ballabio
% Milano Chemometrics and QSAR Research Group
% www.disat.unimib.it/chm

a = param.a;
s = param.s;
m = param.m;
M = param.M;
pret_type = param.pret_type;

if strcmp(pret_type,'cent')
    for j=1:size(X_scal,2)
        X(:,j) = X_scal(:,j) + a(j);
    end
elseif strcmp(pret_type,'scal')
    for j=1:size(X_scal,2)
        X(:,j) = X_scal(:,j)*s(j);
    end
elseif strcmp(pret_type,'auto')
    for j=1:size(X_scal,2)
        X(:,j) = X_scal(:,j)*s(j) + a(j);
    end
elseif strcmp(pret_type,'rang')
    for j=1:size(X_scal,2)
        X(:,j) = X_scal(:,j)*(M(j) - m(j)) + m(j);
    end
else
    X = X_scal;
end
