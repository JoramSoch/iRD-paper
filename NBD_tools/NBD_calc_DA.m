function [DA, BA] = NBD_calc_DA(z, Zp, iZ)
% _
% Calculate Decoding Accuracy (for Classification)
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 10/07/2020, 20:55


% get model dimensions
n = numel(z);
v = size(Zp,2);
q = 2;

% set indices if necessary
if nargin < 3 || isempty(iZ)
    iZ = true(n,1);
end;

% restrict behavioral responses
z  = z(iZ);
Zp = Zp(iZ,:);
tZ = sum(iZ);

% calculate decoding accuracy
Zc = (Zp>(1/q));
DA = (1/tZ) * sum(repmat(z,[1 v])==Zc,1);
BA = (1/2) * mean(repmat(z(z==0),[1 v])==Zc(z==0,:),1) + ...
     (1/2) * mean(repmat(z(z==1),[1 v])==Zc(z==1,:),1);