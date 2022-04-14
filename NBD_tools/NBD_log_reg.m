function [Y_pred, B_est] = NBD_log_reg(Y1, X1, Y2, X2)
% _
% Estimate and Predict using Logistic Regression
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 17/07/2020, 15:06


% get model dimensions
n1 = size(Y1,1);
n2 = size(Y2,1);
p  = size(X1,2);
v  = size(Y1,2);

% estimate parameters
B_est = zeros(p,v);
for j = 1:v
    B_est(:,j) = ME_log_reg_IRLS(Y1(:,j), X1, 0);
end;

% predict data
Y_pred = zeros(n2,v);
for j = 1:v
    Y_pred(:,j) = 1./(1 + exp(-X2*B_est(:,j)));
end;