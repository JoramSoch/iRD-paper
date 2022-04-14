function [Y_pred, B_est] = NBD_lin_reg(Y1, X1, Y2, X2)
% _
% Estimate and Predict using Linear Regression
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 04/12/2019, 09:50


% estimate parameters
B_est = (X1'*X1)^(-1) * X1'*Y1;

% predict data
Y_pred = X2*B_est;