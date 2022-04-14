function NARPS_stats_JS(data_dir, subj_id, MS_name, GLM_names, run)
% _
% Statistical modelling of post-processed NARPS data.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 28/07/2020, 15:49 / 30/07/2020, 13:06


%%% Step 0: Prepare stats parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% experiment parameters
S  = 4;                         % number of sessions
n  = 453;                       % number of scans
TR = 1;                         % repetition time

% create model folders
subj_dir = strcat(data_dir,subj_id,'/');
stat_dir = strcat(subj_dir,'glms/');
if ~exist(stat_dir,'dir'), mkdir(stat_dir); end;
MS_dir   = strcat(stat_dir,'glms-',MS_name,'/');
if ~exist(MS_dir,'dir'), mkdir(MS_dir); end;
GLM_dirs = cell(size(GLM_names));
for k = 1:numel(GLM_names)
    GLM_dirs{k} = strcat(MS_dir,'glm-',GLM_names{k},'/');
    if ~exist(GLM_dirs{k},'dir'), mkdir(GLM_dirs{k}); end;
end;

% create names, onsets, durations
create_onset_files(data_dir, subj_id, MS_name);

% create multiple regressors
create_mult_regs(data_dir, subj_id, [27:32]);


%%% Step 1: Create modelling batches %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find files
func_dir   = strcat(data_dir,'derivatives/fMRIprep/',subj_id,'/func/');
func_files = dir(strcat(func_dir,'w',subj_id,'*space-MNI152NLin2009cAsym_preproc.nii'));
mask_file  = dir(strcat(func_dir,'w',subj_id,'*space-MNI152NLin2009cAsym_GM.nii'));

% create batch: model-independent
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir            = GLM_dirs(1);
matlabbatch{1}.spm.stats.fmri_spec.timing.units   = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
for i = 1:S
    files = cell(n,1);
    for j = 1:n
        files{j} = strcat(func_files(i).folder,'/',func_files(i).name,',',num2str(j));
    end;
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans     = files;
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond      = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi     = {'multi'};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress   = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = {strcat(func_dir,subj_id,sprintf('_task-MGT_run-%02d_bold_confounds.mat',i))};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf       = 128;
end;
matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask             = {strcat(mask_file.folder,'/',mask_file.name,',1')};
matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1)         = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical  = 1;

% create batch: model-dependent
for k = 1:numel(GLM_names)
    matlabbatch{1}.spm.stats.fmri_spec.dir = GLM_dirs(k);
    for i = 1:S
        file{1} = strcat(MS_dir,subj_id,'_glm-',GLM_names{k},'_run-',sprintf('%02d',i),'_onsets.mat');
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = file;
    end;
	filename = strcat(MS_dir,subj_id,'_glm-',GLM_names{k},'_design.mat');
    save(filename,'matlabbatch');
end;


%%% Step 2: Execute modelling batches %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run batches
for k = 1:numel(GLM_names)
    filename = strcat(MS_dir,subj_id,'_glm-',GLM_names{k},'_design.mat');
    if run, spm_jobman('run',filename); end;
end;