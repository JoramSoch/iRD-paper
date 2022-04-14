function [X_pred, DA_in] = NBD_SV_reg_SL(X1, Y1, X2, Y2, SLs)
% _
% Estimate and Predict using Support Vector Regression in Searchlights
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 09/09/2020, 11:34
%  Last edit: 09/09/2020, 11:34


% get model dimensions
p  = size(X1,2);                % number of variables
v  = size(Y1,2);                % number of time series
n1 = size(Y1,1);                % number of training data points
n2 = size(Y2,1);                % number of test data points
d  = floor(v/100);

% create intercepts
o1 = ones(n1,1);                % constant regressor (training)
o2 = ones(n2,1);                % constant regressor (test)

% searchlight-based SVR
X_fitt = zeros(n1,p,v);
X_pred = zeros(n2,p,v);
spm_progress_bar('Init', 100, 'Estimate support vector regressions...', '');
for j = 1:v
    % train and predict
    Y1j = [Y1(:,SLs{j}), o1];
    Y2j = [Y2(:,SLs{j}), o2];
    for k = 1:p
        % calibrate SVM
        svm1 = fitrsvm(Y1j, X1(:,k));
        % training set: in-sample fitting
        X_fitt(:,k,j) = predict(svm1, Y1j);
        % test set: out-of-sample predictions
        X_pred(:,k,j) = predict(svm1, Y2j);
    end;
    % update progress bar
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
X_fitt = permute(X_fitt,[1 3 2]);
X_pred = permute(X_pred,[1 3 2]);
spm_progress_bar('Clear');
clear svm1

% calculate in-sample accuracies
DA_in = zeros(p,v);
for k = 1:p
    DA_in(k,:) = NBD_calc_DP(X1(:,k), X_fitt(:,:,k));
end;
clear X_fitt