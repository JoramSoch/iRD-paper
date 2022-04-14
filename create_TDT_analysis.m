function results_mat = create_TDT_analysis(stat_dir, scan_list, mask_file, type, names, vals, CV)
% _
% Create TDT analysis for classification or regression analysis
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 03/03/2021, 12:25; 30/03/2021, 21:28


% Step 1: create cfg structure
cfg             = decoding_defaults;
cfg.analysis    = 'wholebrain';
cfg.results.dir = stat_dir;

% Step 2: specify files and values
cfg.files.name  = scan_list;
cfg.files.mask  = mask_file;
if strcmp(type,'SVC')
    cfg.files.label              = zeros(size(cfg.files.name));
    cfg.files.labelname          = cell(size(cfg.files.name));
    cfg.files.label(vals==1)     = 1;
    cfg.files.label(vals==2)     = 2;
    cfg.files.labelname(vals==1) = names(1);
    cfg.files.labelname(vals==2) = names(2);
elseif strcmp(type,'SVR')
    cfg.files.label     = vals;
    cfg.files.labelname = repmat(names,[numel(vals),1]);
end;

% Step 3: specify CV design
cfg.files.chunk = zeros(size(cfg.files.name));
for j = 1:size(CV,2)
    cfg.files.chunk(CV(:,j)==2) = j;
end;
cfg.design = make_design_cv(cfg);
cfg.design.unbalanced_data = 'ok';
cfg.plot_design = 0;
cfg.plot_selected_voxels = 0;

% Step 4: set analysis parameters
cfg.verbose           = 1; % normal output
cfg.results.write     = 2; % no brain images
cfg.results.overwrite = 1; % overwrite results
if strcmp(type,'SVC')
    cfg.decoding.method = 'classification_kernel';
    cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 
    cfg.results.output = {'accuracy', 'balanced_accuracy', 'AUC', 'decision_values', 'predicted_labels', 'true_labels'};
elseif strcmp(type,'SVR')
    cfg.decoding.method = 'regression';
    cfg.decoding.train.classification.model_parameters = '-s 4 -t 0 -c 1 -b 0 -q';
    cfg.results.output = {'corr', 'zcorr', 'decision_values', 'predicted_labels', 'true_labels'};
end;

% Step 5: run decoding analysis
[results, cfg] = decoding(cfg);

% return filename
filename = strcat(stat_dir,'results.mat');
save(filename, 'cfg', 'results');
results_mat = filename;