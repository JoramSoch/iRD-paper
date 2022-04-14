function X = NBD_marg2joint(x)
% _
% Transform marginal to joint probabilities assuming independence
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 04/12/2019, 10:50


% get model dimensions
n = size(x,1);                  % number of data points
log2p = size(x,2);              % number of binary factors
p = 2^log2p;                    % number of conditions

% create dummys
D = zeros(log2p,p);
for k = 1:log2p
    D(k,:) = repmat([ones(1,2^(log2p-k)), zeros(1,2^(log2p-k))],[1 2^(k-1)]);
end;

% calculate marginals
M = zeros(n,p,log2p);
for k = 1:log2p
    for i = 1:n
        M(i,:,k) = x(i,k)*(D(k,:)==1) + (1-x(i,k))*(D(k,:)==0);
    end;
end;

% calculate joints
J = prod(M,3);
X = J;
clear J