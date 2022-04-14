function [r, p, MSE, MAE, slo, int] = NBD_calc_DP(z, Zp, iZ)
% _
% Calculate Decoding Precision (for Reconstruction)
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 17/07/2020, 14:52


% get model dimensions
n = numel(z);
v = size(Zp,2);

% set indices if necessary
if nargin < 3 || isempty(iZ)
    iZ = true(n,1);
end;

% restrict behavioral responses
z  = z(iZ);
Zp = Zp(iZ,:);
tZ = sum(iZ);

% calculate decoding accuracy
[r,p] = corr(z, Zp);
MSE   = mean((Zp-repmat(z,[1 v])).^2,1);
MAE   = mean(abs(Zp-repmat(z,[1 v])),1);
slo   = zeros(1,v);
int   = zeros(1,v);
for j = 1:v
    b = polyfit(Zp(:,j), z, 1);
    slo(j) = b(1);
    int(j) = b(2);
end;