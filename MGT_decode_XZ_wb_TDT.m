function MGT_decode_XZ_wb_TDT(data_dir, subj_id, MS_name, GLM_names, ana_id)
% _
% Decode Design X and Behavior Z (TDT)
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 03/03/2021, 12:14 (1/2/3/4d)


%%% Intro: prepare analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n-> Subject "%s", Analysis "%s":\n', subj_id, ana_id);

% check if data already exist
filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','NBD_XZ.mat');
if ~exist(filename,'file')

    % set dimensions
    S   = 4;                        % number of sessions
    tpS = 64;                       % trials per session
    t   = S*tpS;                    % trials in total
    m2s = [mean([5:20]), ...        % means to subtract
           mean([10:40])];

    % load X and Z
    load(strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','X_','glms-',MS_name,'.mat'));
    load(strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','Z_','glms-',MS_name,'.mat'));
    p  = size(X,2);                 % number of design variables
    q  = size(Z,2);                 % number of response variables
    iZ = Z(:,1)~=0;                 % only trials with responses

    % adjust X columns
    if mean(X(:,1)) < 20            % equal range condition
        X(:,1) = X(:,1)-m2s(1);
    else                            % equal indifference condition
        X(:,1) = X(:,1)-m2s(2);
    end;
        X(:,2) = X(:,2)-m2s(1);

    % create X matrices
    Xc = zeros(t,2^p);              % categorical conditions
    Xc(X(:,1)<0 & X(:,2)<0,1) = 1;  % low gain, low loss
    Xc(X(:,1)<0 & X(:,2)>0,2) = 1;  % low gain, high loss
    Xc(X(:,1)>0 & X(:,2)<0,3) = 1;  % high gain, low loss
    Xc(X(:,1)>0 & X(:,2)>0,4) = 1;  % high gain, high loss
    Xp = [X(:,[1 2]), ones(t,1)];   % parametric regressors
    Xr = [ones(t,1), abs(diff(X,1,2))];     % RT predictors

    % create Z matrices
    Zc = zeros(t,2^1);              % categorical responses
    Zc(Z(:,1)==1 | Z(:,1)==2,1) = 1;% (strongly/weakly) reject
    Zc(Z(:,1)==3 | Z(:,1)==4,2) = 1;% (weakly/strongly) accept
    Zp = Z(:,1);                    % parametric responses
    Zr = Z(:,2);                    % reaction times

    % create CV folds
    fprintf('   - CV: Session ');
    is  = [1:t];
    i2  = cell(S,1);           % test set = current session
    i1  = cell(S,1);           % training set = all other sessions
    iZ2 = cell(S,1);           % exclude non-responses
    iZ1 = cell(S,1);           % exclude non-responses
    for g = 1:S
        % test/train indices
        i2{g}  = [(g-1)*tpS+1:g*tpS];
        i1{g}  = setdiff(is, i2{g});
        iZ2{g} = i2{g}( iZ(i2{g}) );
        iZ1{g} = i1{g}( iZ(i1{g}) );
        fprintf('%d, ', g);
    end;
    fprintf('done.\n')

    % X/Z: save matrices
    filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','NBD_XZ.mat');
    save(filename, 'S', 'tpS', 't', 'X', 'Xc', 'Xp', 'Xr', 'Z', 'Zc', 'Zp', 'Zr', 'is', 'i1', 'i2', 'iZ', 'iZ1', 'iZ2');

else
    
    % X/Z: load matrices
    filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','NBD_XZ.mat');
    load(filename);
    
end;


%%% PBM: decode Z from X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if results already exist
filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','NBD_PBM-',ana_id(1),'b','.mat');
if ~exist(filename,'file')

    % PBM: decode Z from X
    fprintf('   - PBM: Session ');
    if strncmp(ana_id,'1',1)
        Zcc = zeros(size(Zc));
        Bcc = zeros(size(Xc,2),size(Zc,2),S);
    end;
    if strncmp(ana_id,'2',1)
        Zpc = zeros(size(Zc,1),1);
        Bpc = zeros(size(Xp,2),S);
    end;
    if strncmp(ana_id,'3',1)
        Zcp = zeros(size(Zp));
        Bcp = zeros(size(Xc,2),S);
    end;
    if strncmp(ana_id,'4',1)
        Zpp = zeros(size(Zp));
        Bpp = zeros(size(Xp,2),S);
    end;
    for g = 1:S
        % Analysis 1a/1b: categorical conditions -> categorical responses
        if strncmp(ana_id,'1',1)
            X1 = Xc(iZ1{g},:);
            X2 = Xc(i2{g},:);
            Z1 = Zc(iZ1{g},:);
            Z2 = Zc(i2{g},:);
            [Zcc(i2{g},:), Bcc(:,:,g)] = NBD_lin_reg(Z1, X1, Z2, X2);
        end;
        % Analysis 2a/2b: parametric regressors -> categorical responses
        if strncmp(ana_id,'2',1)
            X1 = Xp(iZ1{g},:);
            X2 = Xp(i2{g},:);
            Z1 = Zc(iZ1{g},2);
            Z2 = Zc(i2{g},2);
            [Zpc(i2{g},:), Bpc(:,g)] = NBD_log_reg(Z1, X1, Z2, X2);
        end;
        % Analysis 3a/3b: categorical conditions -> parametric responses
        if strncmp(ana_id,'3',1)
            X1 = Xc(iZ1{g},:);
            X2 = Xc(i2{g},:);
            Z1 = Zp(iZ1{g},:);
            Z2 = Zp(i2{g},:);
            [Zcp(i2{g},:), Bcp(:,g)] = NBD_lin_reg(Z1, X1, Z2, X2);
        end;
        % Analysis 4a/4b: parametric regressors -> parametric responses
        if strncmp(ana_id,'4',1)
            X1 = Xp(iZ1{g},:);
            X2 = Xp(i2{g},:);
            Z1 = Zp(iZ1{g},:);
            Z2 = Zp(i2{g},:);
            [Zpp(i2{g},:), Bpp(:,g)] = NBD_lin_reg(Z1, X1, Z2, X2);
        end;
        fprintf('%d, ', g);
    end;
    fprintf('done.\n'),

    % PBM: save results
    filename = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','NBD_PBM-',ana_id,'.mat');
    if strncmp(ana_id,'1',1), save(filename, 'Zcc', 'Bcc'); end;
    if strncmp(ana_id,'2',1), save(filename, 'Zpc', 'Bpc'); end;
    if strncmp(ana_id,'3',1), save(filename, 'Zcp', 'Bcp'); end;
    if strncmp(ana_id,'4',1), save(filename, 'Zpp', 'Bpp'); end;
    
else
    
    % PBM: copy results
    srce_file = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','NBD_PBM-',ana_id(1),'b','.mat');
    dest_file = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','NBD_PBM-',ana_id,'.mat');
    copyfile(srce_file, dest_file, 'f');
    load(dest_file);
    
end;
fprintf('\n');


%%% CRD: decode Z from Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analysis 1d/2d: measured data -> categorical responses (SVC)
if strcmp(ana_id,'1d') || strcmp(ana_id,'2d')
    type  = 'SVC';
    names ={'reject'; 'accept'};
    vals  = (Zc(:,1)==1)*1 + (Zc(:,2)==1)*2;
    CV    = kron(eye(S),ones(tpS,1))+1;
end;

% Analysis 3d/4d: measured data -> parametric responses (SVR)
if strcmp(ana_id,'3d') || strcmp(ana_id,'4d')
    type  = 'SVR';
    names ={'favorability'};
    vals  = Zp;
    CV    = kron(eye(S),ones(tpS,1))+1;
end;

% load parameter estimates
GLM_dir   = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','glm-',GLM_names{1},'/');
load(strcat(GLM_dir,'ITEM_est_1st_lvl','/','GLM1.mat'));

% prepare TDT analysis
stat_dir  = strcat(GLM_dir,'TDT_',type,'_fav','/');
mask_file = strcat(GLM_dir,'mask.nii');
scan_list = cell(t,1);
for i = 1:S
    for j = 1:GLM1.t(i)
        filename = strcat(GLM1.swd, GLM1.Vgamma(GLM1.Sess(i).t(j)).fname);
        scan_list{(i-1)*tpS+j,1} = filename;
    end;
end;

% restrict to actual responses
scan_list = scan_list(iZ);
vals      = vals(iZ);
CV        = CV(iZ,:);

% check if results already exist
fprintf('   - CRD: TDT analysis:\n');
if ~exist(strcat(stat_dir,'results.mat'),'file')
    create_TDT_analysis(stat_dir, scan_list, mask_file, type, names, vals, CV);
end;
fprintf('\n');


%%% NBD: decode X from Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analysis 1d/3d: measured data -> categorical responses (SVC)
if strcmp(ana_id,'1d') || strcmp(ana_id,'2d')
    type  = 'SVC';
    names ={'gain-low',  'loss-low';
            'gain-high', 'loss-high'};
    vals  =[(Xc(:,1)==1 | Xc(:,2)==1)*1 + (Xc(:,3)==1 | Xc(:,4)==1)*2, ...
            (Xc(:,1)==1 | Xc(:,3)==1)*1 + (Xc(:,2)==1 | Xc(:,4)==1)*2];
    CV    = kron(eye(S),ones(tpS,1))+1;
end;

% Analysis 2d/4d: measured data -> parametric responses (SVR)
if strcmp(ana_id,'2d') || strcmp(ana_id,'4d')
    type  = 'SVR';
    names ={'gain', 'loss'};
    vals  =[Xp(:,1), Xp(:,2)];
    CV    = kron(eye(S),ones(tpS,1))+1;
end;

% load parameter estimates
GLM_dir   = strcat(data_dir,subj_id,'/glms/','glms-',MS_name,'/','glm-',GLM_names{2},'/');
load(strcat(GLM_dir,'ITEM_est_1st_lvl','/','GLM1.mat'));

% prepare TDT analysis
stat_dir  ={strcat(GLM_dir,'TDT_',type,'_gain','/'), strcat(GLM_dir,'TDT_',type,'_loss','/')};
mask_file = strcat(GLM_dir,'mask.nii');
scan_list = cell(t,1);
for i = 1:S
    for j = 1:GLM1.t(i)
        filename = strcat(GLM1.swd, GLM1.Vgamma(GLM1.Sess(i).t(j)).fname);
        scan_list{(i-1)*tpS+j,1} = filename;
    end;
end;

% check if results already exist
fprintf('   - NBD: TDT analyses:\n');
if ~exist(strcat(stat_dir{1},'results.mat'),'file')
    create_TDT_analysis(stat_dir{1}, scan_list, mask_file, type, names(:,1), vals(:,1), CV);
end;
if ~exist(strcat(stat_dir{2},'results.mat'),'file')
    create_TDT_analysis(stat_dir{2}, scan_list, mask_file, type, names(:,2), vals(:,2), CV);
end;
fprintf('\n');