% The mixed gambling task (MGT) indirect response decoding (iRD) analysis
% _
% This script performs data analyses underlying the results reported in the
% paper by Soch & Haynes (2022): "Indirect decoding of behavior from brain
% signals via reconstruction of experimental conditions".
% 
% For more info, see:
% - https://github.com/JoramSoch/ITEM/blob/master/README.md
% - https://www.biorxiv.org/content/10.1101/2022.03.24.485588v1


%%% Step 0: prepare data analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run preparatory scripts
addpath(strcat(pwd,'/NBD_tools/'));
project_directories
extract_subjects
extract_trials
NARPS_model_spaces

% load general variables
load project_directories.mat
load extract_subjects.mat
load extract_trials.mat
load MS_files/MS_MGT1.mat
num_subj = numel(subj_ids);
% subj_ids = {'sub-001', ..., 'sub-124'}';


%%% Step 1: preprocessing and statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% perform preproc & stats
for i = 1:num_subj
    fprintf('%s\n', subj_ids{i});
  % NARPS_rename_unzf(data_dir, subj_ids{i});
    NARPS_postproc_JS(data_dir, subj_ids{i}, [1:3], false, true);
    NARPS_stats_JS(data_dir, subj_ids{i}, MS_name, GLM_names, true);
end;


%%% Step 2: data extraction and decoding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract X/Y/Z data
for i = 1:num_subj
    MGT_estimate_Y(data_dir, subj_ids{i}, MS_name, GLM_names);
    MGT_extract_XYZ(data_dir, subj_ids{i}, MS_name, GLM_names);
end;

% perform TDT analyses
ana_ids = {'1d', '2d', '3d', '4d'};
for i = 1:num_subj
    for j = 1:numel(ana_ids)
        MGT_decode_XZ_wb_TDT(data_dir, subj_ids{i}, MS_name, GLM_names, ana_ids{j});
        MGT_show_results_wb_TDT(data_dir, subj_ids{i}, MS_name, GLM_names, ana_ids{j}, 0);
    end;
end;


%%% Step 3: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot figures
Figure_4_5_6
Figure_S3_S4_S5
Figure_S6_S7_S8