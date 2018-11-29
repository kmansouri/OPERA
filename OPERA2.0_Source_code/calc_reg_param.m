function [res,RSS] = calc_reg_param(y_real,y_calc)

% calculation of regression parameters
%
% [res] = reg_eval(y_real,y_calc)
%
% INPUT
% y_real:   vector of experimental responses (samples x 1)
% y_calc:   vector of calculated responses (samples x 1)
% 
% OUTPUT
% res is a structure array with fields:
% R2:       percentage of explained variance
% R:        coefficient of multiple correlation
% RMSEC:    Root Mean Squared Error of Calibration
% 
% version 1.0 - september 2009
% Davide Ballabio
% Milano Chemometrics and QSAR Research Group
% www.disat.unimib.it/chm

n = length(y_real);

% TSS, Total Sum of Squares
TSS = sum((y_real - mean(y_real)).^2);

% RSS, Residual Sum of Squares
RSS = sum((y_real - y_calc).^2);

% R2, percentage of explained variance
R2 = 1 - RSS/TSS;

% R, coefficient of multiple correlation
R = R2^0.5;

% RMSEC, Root Mean Squared Error of Calibration
% also called SDEC, Standard Deviation Error in Calibration
RMSEC = (RSS/n)^0.5;

res.R2 = R2;
res.R = R;
res.RMSEC = RMSEC;