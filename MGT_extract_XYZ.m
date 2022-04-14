function MGT_extract_XYZ(data_dir, subj_id, MS_name, GLM_names)
% _
% Extract Design X, Data Y and Behavior Z
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 22/07/2020, 15:34


% get directory
tool_dir = pwd;

% extract data Y
for j = 1:numel(GLM_names)
    if strcmp(GLM_names{j},'resp') || strcmp(GLM_names{j},'stim')
        % load SPM
        load(strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','glm-',GLM_names{j},'/','SPM.mat'));
        load(strcat(SPM.swd,'/ITEM_est_1st_lvl/','GLM1.mat'));
        % load mask
        [M, m_dim, m_ind] = MA_load_mask(SPM);
        % load gammas
        Y = ITEM_load_gammas(SPM, m_ind);
        V = blkdiag(GLM1.Sess.U);
        % reduce gammas
        i = [];
        for k = 1:numel(GLM1.Sess)
            i = [i, GLM1.Sess(k).t(1:GLM1.t(k))];
        end;
        Y = Y(i,:);
        V = V(i,i);
        % save gammas
        filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/Y_','glm-',GLM_names{j},'.mat');
        save(filename, 'Y', 'V', 'i');
    end;
end;

% load subjects and trials
load(strcat(tool_dir,'/','extract_subjects.mat'));
load(strcat(tool_dir,'/','extract_trials.mat'));
i = find(strcmp(subj_ids,subj_id));

% extract design X
X = vertcat(Subj(i).Sess.T);
X = X(:,4:5);
l = {'gain', 'loss'};
filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','X_','glms-',MS_name,'.mat');
save(filename, 'X', 'l');

% extract behavior Z
Z = vertcat(Subj(i).Sess.T);
Z = Z(:,6:7);
l = {'resp', 'RT'};
filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','Z_','glms-',MS_name,'.mat');
save(filename, 'Z', 'l');

% return to tools
cd(tool_dir)