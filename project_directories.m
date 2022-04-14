% OpenNeuro ds001734: project directories
% _
% This script generates project directories for the OpenNeuro ds001734.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 09/07/2020, 14:45


% set directories
stud_dir = 'I:\joram\NBD\DataSets\OpenNeuro_ds001734\';
data_dir = strcat(stud_dir,'data\');
ownc_dir = 'C:\Joram\ownCloud\BCCN\NBD\DataSets\OpenNeuro_ds001734\';
tool_dir = strcat(ownc_dir,'tools\');

% save directories
save('project_directories.mat','stud_dir','data_dir','tool_dir');