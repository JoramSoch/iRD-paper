function NARPS_rename_unzf(data_dir, subj_id)
% _
% Rename unzipped functional pre-processed NARPS data.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 09/09/2020, 14:02


%%% Step 0: Specify renaming parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimensions of experiment
S = 4;                          % number of sessions
n = 453;                        % number of scans

% subject's functional MRI directory
func_dir = strcat(data_dir,'derivatives/fMRIprep/',subj_id,'/func/');

% unzipped functional MRI scans
unzf_pref =  'vol0000_xform-00000_merged';
unzf_suff = {'', '_1', '_2', '_3'};
unzf_ext  =  '.nii';

% NOTE: These settings assume that the files
%     sub-001_task-MGT_run-{01/02/03/04}_bold_space-MNI152NLin2009cAsym_preproc.nii.gz
% were unpacked using 7-Zip into the files
%     vol0000_xform-00000_merged{/_1/_2/_3}.nii
% and have to be renamed into
%     sub-001_task-MGT_run-01_bold_space-MNI152NLin2009cAsym_preproc.nii
% before post-processing.


%%% Step 1: Rename unzipped functional scans %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find destination files
func_files = dir(strcat(func_dir,subj_id,'*space-MNI152NLin2009cAsym_preproc.nii.gz'));

% move source files
for i = 1:S
    srce_file = strcat(func_dir, unzf_pref, unzf_suff{i}, unzf_ext);
	dest_file = strcat(func_files(i).folder,'/',func_files(i).name(1:end-3));
	movefile(srce_file, dest_file, 'f');
end;