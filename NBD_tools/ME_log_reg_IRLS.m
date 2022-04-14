function beta_est = ME_log_reg_IRLS(y, X, l)
% _
% Estimate Logistic regression using Iterative Reweighted Least Squares
% FORMAT beta_est = ME_log_reg_IRLS(y, X)
% 
%     y - an n x 1 vector of observed categories, coded as 0s and 1s
%     X - an n x p matrix of predictor variables, positive or negative
%     l - a scalar, the regularization parameter for IRLS
% 
%     beta_est - a p x 1 vector of logistic regression coefficients
% 
% FORMAT beta_est = ME_log_reg_IRLS(y, X) returns iterative reweighted
% least sqaures [1] estimates for the coefficients of a logistic
% regression model with design matrix X and data vector y.
% 
% References:
% [1] Bishop CM (2006): "Pattern Recognition and Machine Learning",
%     chs. 4.3.3/4.3.4, pp. 205-208, eqs. 4.89/4.99.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 06/03/2020, 18:58
%  Last edit: 10/03/2020, 09:45


% Get model dimensions
%-------------------------------------------------------------------------%
n = size(y,1);                  % number of data points
p = size(X,2);                  % number of parameters

% Initialize parameters
%-------------------------------------------------------------------------%
b = zeros(p,1);
Ip= eye(p);
conv = 1e-4;
diff = 1;

% Estimate parameters
%-------------------------------------------------------------------------%
while diff > conv
    % calculate class probabilities
    p = 1 ./ (1 + exp(-X*b));
    if any(p==1) || any(p==0), break; end;
    % establish reweighting matrix
    R = diag(p.*(1-p));
    % update regression coefficients
    D = (X'*R*X + l*Ip)^(-1) * X'*(p-y);
    b = b - D;
    diff = norm(D);
end;

% Output parameters
%-------------------------------------------------------------------------%
beta_est = b;