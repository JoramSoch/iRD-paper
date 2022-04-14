% OpenNeuro ds001734: extract subjects
% _
% This script extracts all subjects.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 09/07/2020, 14:43


% load directories
load project_directories.mat

% read subject IDs
[num, txt, raw] = xlsread('subjects.xls');
subjects = raw(2:end,:);
subj_ids = subjects(:,1);
num_subj = numel(subj_ids);

% 1st column: subject index
% 2nd column: subject group (1 = equal range, 2 = equal indifference)
% 3rd column: gender (1 = m, 2 = f)
% 4th column: age [yrs]

% create subjects matrix
S = zeros(num_subj,4);
for i = 1:num_subj
    S(i,1) = str2num(subj_ids{i}(5:end));
end;
S(strcmp(subjects(:,2),'equalRange'),2)        = 1;
S(strcmp(subjects(:,2),'equalIndifference'),2) = 2;
S(strcmp(subjects(:,3),'M'),3) = 1;
S(strcmp(subjects(:,3),'F'),3) = 2;
S(:,4) = cell2mat(subjects(:,4));
clear subjects
Subjects = S;

% save all subjects
fprintf('\n-> Save %d subjects... ', num_subj);
save('extract_subjects.mat', 'Subjects', 'subj_ids');
fprintf('successful!\n\n');