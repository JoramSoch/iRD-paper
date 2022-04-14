function create_mult_regs(data_dir, subj_id, col_inds)
% _
% Create multiple regressors file for MGT data set.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 28/07/2020, 15:53


% set number of sessions
S = 4;

% create multiple regressors
for i = 1:S
	
	% load confounds file
	filename = strcat(data_dir,'derivatives/fmriprep/',subj_id,'/func/',subj_id,sprintf('_task-MGT_run-%02d_bold_confounds.tsv',i));
	[num, hdr, raw] = tsvread(filename);
	if isempty(col_inds)        % all confounds
		R = num(2:end,:);
    else
        try                     % selected confounds
            R = num(2:end,col_inds);
        catch                   % movement parameters
            R = num(2:end,[end-5:end]);
        end;
	end;
		
	% save regressor matrix
	filename = strcat(data_dir,'derivatives/fmriprep/',subj_id,'/func/',subj_id,sprintf('_task-MGT_run-%02d_bold_confounds.mat',i));
	save(filename, 'R');
	
end;