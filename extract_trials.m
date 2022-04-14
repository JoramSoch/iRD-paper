% OpenNeuro ds001734: extract trials
% _
% This script extracts all trials from all sessions of all subjects.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 09/07/2020, 13:52


% set numbers
N = 108;                        % number of subjects in study
S = 4;                          % number of sessions per subject
t = 64;                         % number of trials per session

% load directories
load project_directories.mat

% read subject IDs
[num, txt, raw] = xlsread('subjects.xls');
subjects = raw(2:end,:);
subj_ids = subjects(:,1);
num_subj = numel(subj_ids);
clear subjects

% for all subjects
for i = 1:N
    
    if i == 1, fprintf('\n'); end;
    fprintf('-> Subject %s: Session ', subj_ids{i});
    
    % for all sessions
    for j = 1:S
        
        fprintf('%d, ', j);
        
        % load events file
        filename = strcat(data_dir,subj_ids{i},'/func/',subj_ids{i},'_task-MGT_run-',sprintf('%02d',j),'_events.tsv');
        [num, hdr, raw] = tsvread(filename);
        trials = num(2:end,1:end-1);
        resps  = raw(2:end,end);
        
        % 1st column: trial index
        % 2nd column: trial onset [s]
        % 3rd column: trial duration [s]
        % 4th column: gain value
        % 5th column: loss value
        % 6th column: response (1-4)
        % 7th column: reaction time [s]
        
        % create trials matrix
        T = zeros(t,7);
        T(:,1) = [1:t]';
        T(:,2) = trials(:,1);
        T(:,3) = trials(:,2);
        T(:,4) = trials(:,3);
        T(:,5) = trials(:,4);
        T(:,7) = trials(:,5);
        clear trials
        
        % collect trial responses
        T(strcmp(resps,'NoResp'),6)          = 0;
        T(strcmp(resps,'strongly_reject'),6) = 1;
        T(strcmp(resps,'weakly_reject'),6)   = 2;
        T(strcmp(resps,'weakly_accept'),6)   = 3;
        T(strcmp(resps,'strongly_accept'),6) = 4;
        clear resps
        
        % store trials matrix
        Subj(i).Sess(j).T = T;
        
    end;
    
    fprintf('done.\n');
    
end;

% save all trials
fprintf('\n-> Save %d subjects... ', num_subj);
save('extract_trials.mat', 'Subj');
fprintf('successful!\n\n');