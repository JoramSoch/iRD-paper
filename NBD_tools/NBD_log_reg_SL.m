function [X_pred, DA_in] = NBD_log_reg_SL(x1, Y1, x2, Y2, SLs)
% _
% Estimate and Predict using Logistic Regression in Searchlights
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 11/12/2019, 17:00
%  Last edit: 09/09/2020, 11:37


% get model dimensions
v  = size(Y1,2);                % number of time series
n1 = size(Y1,1);                % number of training data points
n2 = size(Y2,1);                % number of test data points
d  = floor(v/100);

% create intercepts
o1 = ones(n1,1);                % constant regressor (training)
o2 = ones(n2,1);                % constant regressor (test)

% searchlight-based logistic regression
X_fitt = zeros(n1,v);
X_pred = zeros(n2,v);
spm_progress_bar('Init', 100, 'Estimate logistic regressions...', '');
for j = 1:v
    % estimate parameters
    Y1j = [Y1(:,SLs{j}), o1];
    Y2j = [Y2(:,SLs{j}), o2];
    b1  = ME_log_reg_IRLS(x1, Y1j, 0);
    % training set: in-sample fitting
    X_fitt(:,j) = 1./(1 + exp(-Y1j*b1));
    % test set: out-of-sample predictions
    X_pred(:,j) = 1./(1 + exp(-Y2j*b1));
    % update progress bar
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
spm_progress_bar('Clear');
clear Y1j Y2j b1

% calculate in-sample accuracies
DA_in = NBD_calc_DA(x1, X_fitt);
clear X_fitt