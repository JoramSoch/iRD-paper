function create_onset_files(data_dir, subj_id, MS_name)
% _
% Create session-wise onset files for MGT data set.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 14/07/2020, 13:25 / 11/08/2020, 11:12


%%% LOAD ONSETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load subjects and trials
load extract_subjects.mat
load extract_trials.mat

% identify subject
i   = find(strcmp(subj_ids,subj_id));       % subject index
S   = numel(Subj(i).Sess);                  % number of sessions
m2s = [mean([5:20]), mean([10:40])];        % means to subtract



%%% MODEL SPACE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(MS_name,'MGT1')

% set names
names1 = {'stim-LG-LL', 'stim-LG-HL', 'stim-HG-LL', 'stim-HG-HL'};
names2 = {'resp-SR', 'resp-WR', 'resp-WA', 'resp-SA', 'resp-none'};

% save onsets
for j = 1:S
    
    % get and adjust trials
    T = Subj(i).Sess(j).T;
    T(:,4) = T(:,4)-m2s(Subjects(i,2));
    T(:,5) = T(:,5)-m2s(1);
    
    %%% MODEL 1: stimuli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % preallocate
    clear names onsets durations
    onsets    = cell(1,numel(names1));
    durations = cell(1,numel(names1));
    
    % get onsets
    onsets{1} = T(T(:,4)<0 & T(:,5)<0,2);
    onsets{2} = T(T(:,4)<0 & T(:,5)>0,2);
    onsets{3} = T(T(:,4)>0 & T(:,5)<0,2);
    onsets{4} = T(T(:,4)>0 & T(:,5)>0,2);
    
    % get durations
    for k = 1:numel(onsets)
        durations{k} = zeros(numel(onsets{k}),1);
    end;
    
    % save onsets
    names = names1;
    filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/',subj_id,'_','glm-stim','_','run-0',num2str(j),'_','onsets.mat');
    save(filename,'names','onsets','durations');
    
    %%% MODEL 1a: trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get durations
    clear durations
    durations{1} = T(T(:,4)<0 & T(:,5)<0,3);
    durations{2} = T(T(:,4)<0 & T(:,5)>0,3);
    durations{3} = T(T(:,4)>0 & T(:,5)<0,3);
    durations{4} = T(T(:,4)>0 & T(:,5)>0,3);
    
    % save onsets
    filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/',subj_id,'_','glm-trl','_','run-0',num2str(j),'_','onsets.mat');
    save(filename,'names','onsets','durations');
    
    %%% MODEL 2: responses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % preallocate
    clear names onsets durations
    onsets    = cell(1,numel(names2));
    durations = cell(1,numel(names2));
    
    % get onsets
    for k = 1:4
        onsets{k} = T(T(:,6)==k,2) + T(T(:,6)==k,7);
    end;
        onsets{5} = T(T(:,6)==0,2) + T(T(:,6)==0,7);
        
    % get durations
    for k = 1:numel(onsets)
        durations{k} = zeros(numel(onsets{k}),1);
    end;
    
    % save onsets
    names = names2;
    incl  = true(size(names));
    for k = 1:numel(onsets)
        if isempty(onsets{k})
            incl(k) = false;
        end;
    end;
    names     = names(incl);
    onsets    = onsets(incl);
    durations = durations(incl);
    filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/',subj_id,'_','glm-resp','_','run-0',num2str(j),'_','onsets.mat');
    save(filename,'names','onsets','durations');
    
end;
    
end;