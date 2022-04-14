function [X_pred, DA_in] = NBD_lin_reg_SL(X1, Y1, V1, X2, Y2, V2, SLs)
% _
% Estimate and Predict using Linear Regression in Searchlights
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 11/12/2019, 17:00
%  Last edit: 09/09/2020, 11:35


% get model dimensions
p  = size(X1,2);                % number of variables
v  = size(Y1,2);                % number of time series
n1 = size(Y1,1);                % number of training data points
n2 = size(Y2,1);                % number of test data points
d  = floor(v/100);

% create intercepts
o1 = ones(n1,1);                % constant regressor (training)
o2 = ones(n2,1);                % constant regressor (test)

% invert covariance
P1 = inv(V1);                   % precision matrix (training)
W2 = sqrtm(inv(V2));            % whitening matrix (test)

% searchlight-based linear regression
X_fitt = zeros(n1,p,v);
X_pred = zeros(n2,p,v);
spm_progress_bar('Init', 100, 'Estimate linear regressions...', '');
for j = 1:v
    % estimate parameters
    Y1j = [Y1(:,SLs{j}), o1];
    Y2j = [Y2(:,SLs{j}), o2];
    B1  = (Y1j'*P1*Y1j)^(-1) * Y1j'*P1*X1;
    % training set: in-sample fitting
    X_fitt(:,:,j) = Y1j*B1;
    % test set: out-of-sample predictions
    X_pred(:,:,j) = W2*Y2j*B1;
    % update progress bar
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
X_fitt = permute(X_fitt,[1 3 2]);
X_pred = permute(X_pred,[1 3 2]);
spm_progress_bar('Clear');
clear Y1j Y2j B1

% calculate in-sample accuracies
if numel(unique(X1))==2         % decoding accuracy
    X_fitt = double(X_fitt(:,:,1)==max(X_fitt,[],3));
    DA_in  = NBD_calc_DA(X1(:,1), X_fitt);
else                            % decoding precision
    DA_in = zeros(p,v);
    for k = 1:p
        DA_in(k,:) = NBD_calc_DP(X1(:,k), X_fitt(:,:,k));
    end;
end;
clear X_fitt