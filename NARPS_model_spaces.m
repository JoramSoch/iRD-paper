%%% MODEL SPACE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MS_name   =  'MGT1';
GLM_names = {'resp', 'stim', 'trl'}';
if ~exist('MS_files','dir'), mkdir('MS_files'); end;
save('MS_files/MS_MGT1.mat', 'MS_name', 'GLM_names');