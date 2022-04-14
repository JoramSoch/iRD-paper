function MGT_estimate_Y(data_dir, subj_id, MS_name, GLM_names)
% _
% Estimate Trial-Wise Response Amplitudes
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 30/07/2020, 14:21


% get directory
tool_dir = pwd;

% estimate data Y
for j = 1:numel(GLM_names)
    if strcmp(GLM_names{j},'resp') || strcmp(GLM_names{j},'stim')
        % calculate trial-wise parameter estimates
        SPM_mat = strcat(data_dir,'/',subj_id,'/glms/','glms-',MS_name,'/','glm-',GLM_names{j},'/','SPM.mat');
        load(SPM_mat);
        ITEM_est_1st_lvl(SPM);
    end;
end;

% return to tools
cd(tool_dir);