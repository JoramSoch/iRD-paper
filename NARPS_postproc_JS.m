function NARPS_postproc_JS(data_dir, subj_id, steps, unzf, run)
% _
% Post-processing of (fMRIprep-based) pre-processed NARPS data.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 28/07/2020, 14:54 / 30/07/2020, 12:21 / 05/08/2020, 15:08


%%% Step 0: Specify post-proc parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimensions of experiment
S = 4;                          % number of sessions
n = 453;                        % number of scans

% pre-processed functional MRI directories
anat_dir = strcat(data_dir,'derivatives/fMRIprep/',subj_id,'/anat/');
func_dir = strcat(data_dir,'derivatives/fMRIprep/',subj_id,'/func/');


%%% Step 1: Unzip pre-processed data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(1,steps)

% find files
anat_file   = dir(strcat(anat_dir,subj_id,'*space-MNI152NLin2009cAsym_preproc.nii.gz'));
func_files1 = dir(strcat(func_dir,subj_id,'*space-MNI152NLin2009cAsym_brainmask.nii.gz'));
func_files2 = dir(strcat(func_dir,subj_id,'*space-MNI152NLin2009cAsym_preproc.nii.gz'));

% create batch
clear matlabbatch
if unzf
    all_files = [anat_file; func_files1; func_files2];
else
    all_files = [anat_file; func_files1];
end;
files = cell(numel(all_files),1);
for j = 1:numel(all_files)
    files{j} = strcat(all_files(j).folder,'/',all_files(j).name);
end;
matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.files = files;

% save batch
filename = strcat(data_dir,subj_id,'/',subj_id,'_postproc_01_gunzip.mat');
save(filename,'matlabbatch');
if run, spm_jobman('run',filename); end;

end;


%%% Step 2: Normalize pre-processed data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(2,steps)

% find files
anat_file   = dir(strcat(anat_dir,subj_id,'*space-MNI152NLin2009cAsym_preproc.nii'));
func_files1 = dir(strcat(func_dir,subj_id,'*space-MNI152NLin2009cAsym_brainmask.nii'));
func_files2 = dir(strcat(func_dir,subj_id,'*space-MNI152NLin2009cAsym_preproc.nii'));

% create batch: files
clear matlabbatch
file  = cell(1,1);
files = cell(S+S*n,1);
file{1} = strcat(anat_file.folder,'/',anat_file.name,',1');
for i = 1:S
    files{0+i} = strcat(func_files1(i).folder,'/',func_files1(i).name,',1');
end;
for i = 1:S
    for j = 1:n
        files{S+(i-1)*n+j} = strcat(func_files2(i).folder,'/',func_files2(i).name,',',num2str(j));
    end;
end;
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol      = file;
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = files;

% create batch: parameters
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg  = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm      = {strcat(spm('dir'),'/tpm/','TPM.nii')};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg   = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg      = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm     = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp     = 3;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb       = [-78, -112, -70; 78, 76, 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox      = [3 3 3];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp   = 4;

% save batch
filename = strcat(data_dir,subj_id,'/',subj_id,'_postproc_02_normalize.mat');
save(filename,'matlabbatch');
if run, spm_jobman('run',filename); end;

end;


%%% Step 3: Calculate grey matter mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(3,steps)

% find files
func_files = dir(strcat(func_dir,'w',subj_id,'*space-MNI152NLin2009cAsym_brainmask.nii'));

% create batch: files
clear matlabbatch
files = cell(S+1,1);
for i = 1:S
    files{i} = strcat(func_files(i).folder,'/',func_files(i).name,',1');
end;
files{S+1}   = strcat(spm('dir'),'/tpm/','TPM.nii',',1');
matlabbatch{1}.spm.util.imcalc.input  = files;
matlabbatch{1}.spm.util.imcalc.output = strcat('w',subj_id,'_task-MGT_runs-01234_bold_space-MNI152NLin2009cAsym_GM.nii');
matlabbatch{1}.spm.util.imcalc.outdir = {func_dir};

% create batch: parameters
str = '';
for i = 1:S, str = strcat(str,'(i',num2str(i),'~=0).*'); end;
str = strcat(str,'(i5>1/3)');
matlabbatch{1}.spm.util.imcalc.expression     = str;
matlabbatch{1}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;

% save batch
filename = strcat(data_dir,subj_id,'/',subj_id,'_postproc_03_GM_mask.mat');
save(filename,'matlabbatch');
if run, spm_jobman('run',filename); end;

end;