% OpenNeuro ds001734: behavioral models
% _
% This script fits behavioral models to subjects' responses:
% - Xc -> Zc: multivariate regression
% - Xp -> Zc: logistic regression
% - Xc -> Zp: linear regression
% - Xp -> Zp: linear regression
% - Xc -> Zr: linear regression
% - Xp -> Zr: linear regression
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 10/07/2020, 16:15 (V1) / 17/07/2020, 15:59 (V2) /
%         31/07/2020, 15:12 (V3) / 13/08/2020, 15:16 (V4)


clear
close all

%%% Step 0: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load behavioral data
load extract_subjects.mat
load extract_trials.mat

% prepare data analysis
N = numel(Subj);                % number of subjects
S = numel(Subj(1).Sess);        % number of sessions
t = size(Subj(1).Sess(1).T,1);  % number of trials
p = 4;
q = 2;

% specify analysis parameters
m2s = [mean([5:20]), mean([10:40])];        % means to subtract
grl = {'equal range','equal indifference'}; % group labels
col =  'rgb';                               % plotting colors
par = [2, 1; 3, 1; 3, 2];                   % parameter indices


%%% Step 1: extract matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for all subjects
for i = 1:N
    
    % for all sessions
    for j = 1:S
        
        % extract trial matrix
        T = Subj(i).Sess(j).T;
        X = zeros(t,p);
        Z = zeros(t,q);
        
        % adjust gains and losses
        Xr = [ones(t,1), abs(diff(T(:,[4 5]),1,2))];% RT predictors
        if mean(T(:,4)) < 20                % equal range condition
            T(:,4) = T(:,4)-m2s(1);
        else                                % equal indifference condition
            T(:,4) = T(:,4)-m2s(2);
        end;
            T(:,5) = T(:,5)-m2s(1);
        
        % create design matrix
        X(T(:,4)<0 & T(:,5)<0,1) = 1;       % low gain, low loss
        X(T(:,4)<0 & T(:,5)>0,2) = 1;       % low gain, high loss
        X(T(:,4)>0 & T(:,5)<0,3) = 1;       % high gain, low loss
        X(T(:,4)>0 & T(:,5)>0,4) = 1;       % high gain, high loss
        Xp = [T(:,4), T(:,5), ones(t,1)];   % parametric regressors
        
        % create response matrix
        Z(T(:,6)==1 | T(:,6) == 2,1) = 1;   % (strongly/weakly) reject
        Z(T(:,6)==3 | T(:,6) == 4,2) = 1;   % (weakly/strongly) accept
        Zp = T(:,6);                        % parametric responses
        Zr = T(:,7);                        % reaction times
        
        % store matrices
        Subj(i).Sess(j).Xc = X;
        Subj(i).Sess(j).Zc = Z;
        Subj(i).Sess(j).Xp = Xp;
        Subj(i).Sess(j).Zp = Zp;
        Subj(i).Sess(j).Xr = Xr;
        Subj(i).Sess(j).Zr = Zr;
        
    end;
    
    % concatenate matrices
    Subj(i).Xc = vertcat(Subj(i).Sess.Xc);
    Subj(i).Zc = vertcat(Subj(i).Sess.Zc);
    Subj(i).Xp = vertcat(Subj(i).Sess.Xp);
    Subj(i).Zp = vertcat(Subj(i).Sess.Zp);
    Subj(i).Xr = vertcat(Subj(i).Sess.Xr);
    Subj(i).Zr = vertcat(Subj(i).Sess.Zr);
    
end;


%%% Step 2: define CV folds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for all subjects
for i = 1:N
    
    % get trials with responses
    iZ = sum(Subj(i).Zc,2)~=0;              % only trials with responses
    is = [1:S*t];                           % all trials across sessions
    
    % for all CV folds
    for j = 1:S
        
        % get training and test indices
        i2  = [(j-1)*t+1:j*t];              % test set = current session
        i1  = setdiff(is,i2);               % training set = other sessions
        iZ2 = i2( iZ(i2) );                 % test trials with responses
        iZ1 = i1( iZ(i1) );                 % training trials with responses
        
        % store training and test indices
        Subj(i).Sess(j).i2  = i2;
        Subj(i).Sess(j).i1  = i1;
        Subj(i).Sess(j).iZ2 = iZ2;
        Subj(i).Sess(j).iZ1 = iZ1;
        
    end;
    
    % store response indices
    Subj(i).iZ = iZ;
    
end;


%%% Step 3: estimate models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for all subjects
for i = 1:N
    
    % preallocate predictions
    Subj(i).Zcc = zeros(S*t,2);
    Subj(i).Zpc = zeros(S*t,1);
    Subj(i).Zcp = zeros(S*t,1);
    Subj(i).Zpp = zeros(S*t,1);
    Subj(i).Zcr = zeros(S*t,1);
    Subj(i).Zpr = zeros(S*t,1);
    
    % for all CV folds
    for j = 1:S
        
        % get training and test
        Xc1 = Subj(i).Xc(Subj(i).Sess(j).iZ1,:);
        Xc2 = Subj(i).Xc(Subj(i).Sess(j).iZ2,:);
        Zc1 = Subj(i).Zc(Subj(i).Sess(j).iZ1,:);
        Zc2 = Subj(i).Zc(Subj(i).Sess(j).iZ2,:);
        Xp1 = Subj(i).Xp(Subj(i).Sess(j).iZ1,:);
        Xp2 = Subj(i).Xp(Subj(i).Sess(j).iZ2,:);
        Zp1 = Subj(i).Zp(Subj(i).Sess(j).iZ1,:);
        Zp2 = Subj(i).Zp(Subj(i).Sess(j).iZ2,:);
        Xr1 = Subj(i).Xr(Subj(i).Sess(j).iZ1,:);
        Xr2 = Subj(i).Xr(Subj(i).Sess(j).iZ2,:);
        Zr1 = Subj(i).Zr(Subj(i).Sess(j).iZ1,:);
        Zr2 = Subj(i).Zr(Subj(i).Sess(j).iZ2,:);
        
        % Analysis 1: experimental category -> response category
        [Zcc, Bcc] = NBD_lin_reg(Zc1, Xc1, Zc2, Xc2);
        
        % Analysis 2: parametric regressors -> response category
        [Zpc, Bpc] = NBD_log_reg(Zc1(:,2), Xp1, Zc2(:,2), Xp2);
        
        % Analysis 3: experimental category -> parametric response
        [Zcp, Bcp] = NBD_lin_reg(Zp1, Xc1, Zp2, Xc2);
        
        % Analysis 4: parametric regressors -> parametric response
        [Zpp, Bpp] = NBD_lin_reg(Zp1, Xp1, Zp2, Xp2);
        
        % Analysis 5: experimental category -> reaction time
        [Zcr, Bcr] = NBD_lin_reg(Zr1, Xc1, Zr2, Xc2);
        
        % Analysis 6: parametric regressors -> reaction time
        [Zpr, Bpr] = NBD_lin_reg(Zr1, Xr1, Zr2, Xr2);

        % store results
        Subj(i).Sess(j).Bcc = Bcc;
        Subj(i).Sess(j).Bpc = Bpc;
        Subj(i).Sess(j).Bcp = Bcp;
        Subj(i).Sess(j).Bpp = Bpp;
        Subj(i).Sess(j).Bcr = Bcr;
        Subj(i).Sess(j).Bpr = Bpr;
        Subj(i).Zcc(Subj(i).Sess(j).iZ2,:) = Zcc;
        Subj(i).Zpc(Subj(i).Sess(j).iZ2)   = Zpc;
        Subj(i).Zcp(Subj(i).Sess(j).iZ2)   = Zcp;
        Subj(i).Zpp(Subj(i).Sess(j).iZ2)   = Zpp;
        Subj(i).Zcr(Subj(i).Sess(j).iZ2)   = Zcr;
        Subj(i).Zpr(Subj(i).Sess(j).iZ2)   = Zpr;
        
    end;
    
    % average patterns
    Subj(i).Bcc = mean(cat(3,Subj(i).Sess.Bcc),3);
    Subj(i).Bpc = mean(horzcat(Subj(i).Sess.Bpc),2);
    Subj(i).Bcp = mean(horzcat(Subj(i).Sess.Bcp),2);
    Subj(i).Bpp = mean(horzcat(Subj(i).Sess.Bpp),2);
    Subj(i).Bcr = mean(horzcat(Subj(i).Sess.Bcr),2);
    Subj(i).Bpr = mean(horzcat(Subj(i).Sess.Bpr),2);
    [Subj(i).DAcc, Subj(i).BAcc] = NBD_calc_DA(Subj(i).Zc(:,2), Subj(i).Zcc(:,2), Subj(i).iZ);
    [Subj(i).DApc, Subj(i).BApc] = NBD_calc_DA(Subj(i).Zc(:,2), Subj(i).Zpc, Subj(i).iZ);
    [Subj(i).rcp, p, Subj(i).MSEcp, MAE] = NBD_calc_DP(Subj(i).Zp, Subj(i).Zcp, Subj(i).iZ);
    [Subj(i).rpp, p, Subj(i).MSEpp, MAE] = NBD_calc_DP(Subj(i).Zp, Subj(i).Zpp, Subj(i).iZ);
    [Subj(i).rcr, p, Subj(i).MSEcr, MAE] = NBD_calc_DP(Subj(i).Zr, Subj(i).Zcr, Subj(i).iZ);
    [Subj(i).rpr, p, Subj(i).MSEpr, MAE] = NBD_calc_DP(Subj(i).Zr, Subj(i).Zpr, Subj(i).iZ);
    
end;

% get decoding accuracies
DAcc = vertcat(Subj.DAcc);
DAcc_mean = [mean(DAcc(Subjects(:,2)==1)), mean(DAcc(Subjects(:,2)==2))];
DAcc_std  = [std(DAcc(Subjects(:,2)==1)),  std(DAcc(Subjects(:,2)==2))];
DApc = vertcat(Subj.DApc);
DApc_mean = [mean(DApc(Subjects(:,2)==1)), mean(DApc(Subjects(:,2)==2))];
DApc_std  = [std(DApc(Subjects(:,2)==1)),  std(DApc(Subjects(:,2)==2))];

% get balanced accuracies
BAcc = vertcat(Subj.BAcc);
BAcc_mean = [mean(BAcc(Subjects(:,2)==1)), mean(BAcc(Subjects(:,2)==2))];
BAcc_std  = [std(BAcc(Subjects(:,2)==1)),  std(BAcc(Subjects(:,2)==2))];
BApc = vertcat(Subj.BApc);
BApc_mean = [mean(BApc(Subjects(:,2)==1)), mean(BApc(Subjects(:,2)==2))];
BApc_std  = [std(BApc(Subjects(:,2)==1)),  std(BApc(Subjects(:,2)==2))];

% get correlations
rcp = vertcat(Subj.rcp);
rcp_mean = [mean(rcp(Subjects(:,2)==1)), mean(rcp(Subjects(:,2)==2))];
rcp_std  = [std(rcp(Subjects(:,2)==1)),  std(rcp(Subjects(:,2)==2))];
rpp = vertcat(Subj.rpp);
rpp_mean = [mean(rpp(Subjects(:,2)==1)), mean(rpp(Subjects(:,2)==2))];
rpp_std  = [std(rpp(Subjects(:,2)==1)),  std(rpp(Subjects(:,2)==2))];
rcr = vertcat(Subj.rcr);
rcr_mean = [mean(rcr(Subjects(:,2)==1)), mean(rcr(Subjects(:,2)==2))];
rcr_std  = [std(rcr(Subjects(:,2)==1)),  std(rcr(Subjects(:,2)==2))];
rpr = vertcat(Subj.rpr);
rpr_mean = [mean(rpr(Subjects(:,2)==1)), mean(rpr(Subjects(:,2)==2))];
rpr_std  = [std(rpr(Subjects(:,2)==1)),  std(rpr(Subjects(:,2)==2))];

% get mean square errors
MSEcp = vertcat(Subj.MSEcp);
MSEcp_mean = [mean(MSEcp(Subjects(:,2)==1)), mean(MSEcp(Subjects(:,2)==2))];
MSEcp_std  = [std(MSEcp(Subjects(:,2)==1)),  std(MSEcp(Subjects(:,2)==2))];
MSEpp = vertcat(Subj.MSEpp);
MSEpp_mean = [mean(MSEpp(Subjects(:,2)==1)), mean(MSEpp(Subjects(:,2)==2))];
MSEpp_std  = [std(MSEpp(Subjects(:,2)==1)),  std(MSEpp(Subjects(:,2)==2))];
MSEcr = vertcat(Subj.MSEcr);
MSEcr_mean = [mean(MSEcr(Subjects(:,2)==1)), mean(MSEcr(Subjects(:,2)==2))];
MSEcr_std  = [std(MSEcr(Subjects(:,2)==1)),  std(MSEcr(Subjects(:,2)==2))];
MSEpr = vertcat(Subj.MSEpr);
MSEpr_mean = [mean(MSEpr(Subjects(:,2)==1)), mean(MSEpr(Subjects(:,2)==2))];
MSEpr_std  = [std(MSEpr(Subjects(:,2)==1)),  std(MSEpr(Subjects(:,2)==2))];

% get average patterns (categorical)
Bcc = cat(3,Subj.Bcc);
Bcc_mean(:,:,1) = mean(Bcc(:,:,Subjects(:,2)==1),3);
Bcc_mean(:,:,2) = mean(Bcc(:,:,Subjects(:,2)==2),3);
Bcc_std(:,:,1)  = std(Bcc(:,:,Subjects(:,2)==1),[],3);
Bcc_std(:,:,2)  = std(Bcc(:,:,Subjects(:,2)==2),[],3);
Bpc = horzcat(Subj.Bpc);

% get average patterns (parametric)
Bcp = horzcat(Subj.Bcp);
Bcp_mean(:,1) = mean(Bcp(:,Subjects(:,2)==1),2);
Bcp_mean(:,2) = mean(Bcp(:,Subjects(:,2)==2),2);
Bcp_std(:,1)  = std(Bcp(:,Subjects(:,2)==1),[],2);
Bcp_std(:,2)  = std(Bcp(:,Subjects(:,2)==2),[],2);
Bpp = horzcat(Subj.Bpp);
Bpp_mean(:,1) = mean(Bpp(:,Subjects(:,2)==1),2);
Bpp_mean(:,2) = mean(Bpp(:,Subjects(:,2)==2),2);
Bpp_std(:,1)  = std(Bpp(:,Subjects(:,2)==1),[],2);
Bpp_std(:,2)  = std(Bpp(:,Subjects(:,2)==2),[],2);

% get average patterns (reactions)
Bcr = horzcat(Subj.Bcr);
Bcr_mean(:,1) = mean(Bcr(:,Subjects(:,2)==1),2);
Bcr_mean(:,2) = mean(Bcr(:,Subjects(:,2)==2),2);
Bcr_std(:,1)  = std(Bcr(:,Subjects(:,2)==1),[],2);
Bcr_std(:,2)  = std(Bcr(:,Subjects(:,2)==2),[],2);
Bpr = horzcat(Subj.Bpr);
Bpr_mean(:,1) = mean(Bpr(:,Subjects(:,2)==1),2);
Bpr_mean(:,2) = mean(Bpr(:,Subjects(:,2)==2),2);
Bpr_std(:,1)  = std(Bpr(:,Subjects(:,2)==1),[],2);
Bpr_std(:,2)  = std(Bpr(:,Subjects(:,2)==2),[],2);

% save all results
save('behavioral_models.mat', 'Subj', 'DAcc', 'DApc', 'BAcc', 'BApc', ...
     'rcp', 'rpp', 'rcr', 'rpr', 'MSEcp', 'MSEpp', 'MSEcr', 'MSEpr', ...
     'Bcc', 'Bpc', 'Bcp', 'Bpp', 'Bcr', 'Bpr');


%%% Step 4: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'MGT: behavioral models (1a)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

subplot(2,4,[1:3]);
bar(DAcc, 'r');
axis([(1-1), (N+1), -0.01, +1.01]);
xlabel('subject index', 'FontSize', 12);
ylabel('decoding accuracy (all sessions)', 'FontSize', 12);
title('multivariate regression: Xc -> Zc', 'FontSize', 16);

subplot(2,4,4); hold on;
bar(DAcc_mean, 'r');
errorbar([1:numel(grl)], DAcc_mean, DAcc_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), -0.01, +1.01]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('decoding accuracy (mean & SD)', 'FontSize', 12);

subplot(2,4,[5:7]);
bar(BAcc, 'g');
axis([(1-1), (N+1), -0.01, +1.01]);
xlabel('subject index', 'FontSize', 12);
ylabel('balanced accuracy (all sessions)', 'FontSize', 12);
title('multivariate regression: Xc -> Zc', 'FontSize', 16);

subplot(2,4,8); hold on;
bar(BAcc_mean, 'g');
errorbar([1:numel(grl)], BAcc_mean, BAcc_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), -0.01, +1.01]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('balanced accuracy (mean & SD)', 'FontSize', 12);

figure('Name', 'MGT: behavioral models (1b)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for h = 1:numel(grl)
    subplot(1,4,(h-1)*2+1);
    imagesc(Bcc_mean(:,:,h));
    caxis([0 1]);
    colorbar;
    set(gca,'XTick',[1:q],'YTick',[1:p]);
    xlabel('response', 'FontSize', 12);
    ylabel('condition', 'FontSize', 12);
    title(sprintf('%s: mean', grl{h}), 'FontSize', 16);
    
    subplot(1,4,(h-1)*2+2);
    imagesc(Bcc_std(:,:,h));
    caxis([0 max(Bcc_std(:))]);
    colorbar;
    set(gca,'XTick',[1:q],'YTick',[1:p]);
    xlabel('response', 'FontSize', 12);
    ylabel('condition', 'FontSize', 12);
    title('SD', 'FontSize', 16);
end;

figure('Name', 'MGT: behavioral models (2a)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

subplot(2,4,[1:3]);
bar(DApc, 'r');
axis([(1-1), (N+1), -0.01, +1.01]);
xlabel('subject index', 'FontSize', 12);
ylabel('decoding accuracy (all sessions)', 'FontSize', 12);
title('logistic regression: Xp -> Zc', 'FontSize', 16);

subplot(2,4,4); hold on;
bar(DApc_mean, 'r');
errorbar([1:numel(grl)], DApc_mean, DApc_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), -0.01, +1.01]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('decoding accuracy (mean & SD)', 'FontSize', 12);

subplot(2,4,[5:7]);
bar(BApc, 'g');
axis([(1-1), (N+1), -0.01, +1.01]);
xlabel('subject index', 'FontSize', 12);
ylabel('balanced accuracy (all sessions)', 'FontSize', 12);
title('logistic regression: Xp -> Zc', 'FontSize', 16);

subplot(2,4,8); hold on;
bar(BApc_mean, 'g');
errorbar([1:numel(grl)], BApc_mean, BApc_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), -0.01, +1.01]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('balanced accuracy (mean & SD)', 'FontSize', 12);

figure('Name', 'MGT: behavioral models (2b)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for l = 1:size(par,1)
    subplot(1,3,l); hold on;
    plot(Bpc(par(l,1),Subjects(:,2)==1), Bpc(par(l,2),Subjects(:,2)==1), '.r');
    plot(Bpc(par(l,1),Subjects(:,2)==2), Bpc(par(l,2),Subjects(:,2)==2), '.b');
    xlim([min(Bpc(par(l,1),:))-1/10*range(Bpc(par(l,1),:)), max(Bpc(par(l,1),:))+1/10*range(Bpc(par(l,1),:))]);
    ylim([min(Bpc(par(l,2),:))-1/10*range(Bpc(par(l,2),:)), max(Bpc(par(l,2),:))+1/10*range(Bpc(par(l,2),:))]);
    axis square;
    set(gca,'Box','On');
    legend(grl, 'Location', 'NorthEast');
    xlabel(['\beta_',num2str(par(l,1))], 'FontSize', 12);
    ylabel(['\beta_',num2str(par(l,2))], 'FontSize', 12);
    if l == 1, title('logistic regression: Xp -> Zc', 'FontSize', 16); end;
end;

figure('Name', 'MGT: behavioral models (3a)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

subplot(2,4,[1:3]);
bar(rcp, 'r');
axis([(1-1), (N+1), -0.01, +1.01]);
xlabel('subject index', 'FontSize', 12);
ylabel('correlation coefficient (all sessions)', 'FontSize', 12);
title('linear regression: Xc -> Zp', 'FontSize', 16);

subplot(2,4,4); hold on;
bar(rcp_mean, 'r');
errorbar([1:numel(grl)], rcp_mean, rcp_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), -0.01, +1.01]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('correlation coefficient (mean & SD)', 'FontSize', 12);

subplot(2,4,[5:7]);
bar(MSEcp, 'g');
axis([(1-1), (N+1), 0, max([max(MSEcp), max(MSEpp)])]);
xlabel('subject index', 'FontSize', 12);
ylabel('mean squared error (all sessions)', 'FontSize', 12);

subplot(2,4,8); hold on;
bar(MSEcp_mean, 'g');
errorbar([1:numel(grl)], MSEcp_mean, MSEcp_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), 0, max([max(MSEcp), max(MSEpp)])]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('mean squared error (mean & SD)', 'FontSize', 12);

figure('Name', 'MGT: behavioral models (3b)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for h = 1:numel(grl)
    subplot(1,2,h); hold on;
    bar([1:size(Bcp,1)], Bcp_mean(:,h), 'b');
    errorbar([1:size(Bcp,1)], Bcp_mean(:,h), Bcp_std(:,h), '.k', 'LineWidth', 2, 'CapSize', 15);
    xlim([(1-1), (size(Bcp,1)+1)]);
    ylim([0, (11/10)*max(max(Bcp_mean+Bcp_std))]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:size(Bcp,1)]);
    xlabel('\beta', 'FontSize', 12);
    ylabel('regression coefficients (mean & SD)', 'FontSize', 12);
    title(grl{h}, 'FontSize', 16);
end;

figure('Name', 'MGT: behavioral models (4a)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

subplot(2,4,[1:3]);
bar(rpp, 'r');
axis([(1-1), (N+1), -0.01, +1.01]);
xlabel('subject index', 'FontSize', 12);
ylabel('correlation coefficient (all sessions)', 'FontSize', 12);
title('linear regression: Xp -> Zp', 'FontSize', 16);

subplot(2,4,4); hold on;
bar(rpp_mean, 'r');
errorbar([1:numel(grl)], rpp_mean, rpp_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), -0.01, +1.01]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('correlation coefficient (mean & SD)', 'FontSize', 12);

subplot(2,4,[5:7]);
bar(MSEpp, 'g');
axis([(1-1), (N+1), 0, max([max(MSEcp), max(MSEpp)])]);
xlabel('subject index', 'FontSize', 12);
ylabel('mean squared error (all sessions)', 'FontSize', 12);

subplot(2,4,8); hold on;
bar(MSEpp_mean, 'g');
errorbar([1:numel(grl)], MSEpp_mean, MSEpp_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), 0, max([max(MSEcp), max(MSEpp)])]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('mean squared error (mean & SD)', 'FontSize', 12);

figure('Name', 'MGT: behavioral models (4b)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for l = 1:size(par,1)
    subplot(1,3,l); hold on;
    plot(Bpp(par(l,1),Subjects(:,2)==1), Bpp(par(l,2),Subjects(:,2)==1), '.r');
    plot(Bpp(par(l,1),Subjects(:,2)==2), Bpp(par(l,2),Subjects(:,2)==2), '.b');
    xlim([min(Bpp(par(l,1),:))-1/10*range(Bpp(par(l,1),:)), max(Bpp(par(l,1),:))+1/10*range(Bpp(par(l,1),:))]);
    ylim([min(Bpp(par(l,2),:))-1/10*range(Bpp(par(l,2),:)), max(Bpp(par(l,2),:))+1/10*range(Bpp(par(l,2),:))]);
    axis square;
    set(gca,'Box','On');
    legend(grl, 'Location', 'NorthEast');
    xlabel(['\beta_',num2str(par(l,1))], 'FontSize', 12);
    ylabel(['\beta_',num2str(par(l,2))], 'FontSize', 12);
    if l == 1, title('linear regression: Xp -> Zp', 'FontSize', 16); end;
end;

figure('Name', 'MGT: behavioral models (5a)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

subplot(2,4,[1:3]);
bar(rcr, 'r');
axis([(1-1), (N+1), -(11/10)*max(abs(rcr)), +(11/10)*max(abs(rcr))]);
xlabel('subject index', 'FontSize', 12);
ylabel('correlation coefficient (all sessions)', 'FontSize', 12);
title('linear regression: Xc -> Zr', 'FontSize', 16);

subplot(2,4,4); hold on;
bar(rcr_mean, 'r');
errorbar([1:numel(grl)], rcr_mean, rcr_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), -(11/10)*max(abs(rcr)), +(11/10)*max(abs(rcr))]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('correlation coefficient (mean & SD)', 'FontSize', 12);

subplot(2,4,[5:7]);
bar(MSEcr, 'g');
axis([(1-1), (N+1), 0, max([max(MSEcr), max(MSEpp)])]);
xlabel('subject index', 'FontSize', 12);
ylabel('mean squared error (all sessions)', 'FontSize', 12);

subplot(2,4,8); hold on;
bar(MSEcr_mean, 'g');
errorbar([1:numel(grl)], MSEcr_mean, MSEcr_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), 0, max([max(MSEcr), max(MSEpp)])]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('mean squared error (mean & SD)', 'FontSize', 12);

figure('Name', 'MGT: behavioral models (5b)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for h = 1:numel(grl)
    subplot(1,2,h); hold on;
    bar([1:size(Bcr,1)], Bcr_mean(:,h), 'b');
    errorbar([1:size(Bcr,1)], Bcr_mean(:,h), Bcr_std(:,h), '.k', 'LineWidth', 2, 'CapSize', 15);
    xlim([(1-1), (size(Bcr,1)+1)]);
    ylim([0, (11/10)*max(max(Bcr_mean+Bcr_std))]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:size(Bcr,1)]);
    xlabel('\beta', 'FontSize', 12);
    ylabel('regression coefficients (mean & SD)', 'FontSize', 12);
    title(grl{h}, 'FontSize', 16);
end;

figure('Name', 'MGT: behavioral models (6a)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

subplot(2,4,[1:3]);
bar(rpr, 'r');
axis([(1-1), (N+1), -(11/10)*max(abs(rpr)), +(11/10)*max(abs(rpr))]);
xlabel('subject index', 'FontSize', 12);
ylabel('correlation coefficient (all sessions)', 'FontSize', 12);
title('linear regression: Xp -> Zr', 'FontSize', 16);

subplot(2,4,4); hold on;
bar(rpr_mean, 'r');
errorbar([1:numel(grl)], rpr_mean, rpr_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), -(11/10)*max(abs(rpr)), +(11/10)*max(abs(rpr))]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('correlation coefficient (mean & SD)', 'FontSize', 12);

subplot(2,4,[5:7]);
bar(MSEpr, 'g');
axis([(1-1), (N+1), 0, max([max(MSEcp), max(MSEpr)])]);
xlabel('subject index', 'FontSize', 12);
ylabel('mean squared error (all sessions)', 'FontSize', 12);

subplot(2,4,8); hold on;
bar(MSEpr_mean, 'g');
errorbar([1:numel(grl)], MSEpr_mean, MSEpr_std, '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (2.5), 0, max([max(MSEcp), max(MSEpr)])]);
set(gca,'Box','On');
set(gca,'XTick',[1:2],'XTickLabel',grl);
xlabel('subject group', 'FontSize', 12);
ylabel('mean squared error (mean & SD)', 'FontSize', 12);

figure('Name', 'MGT: behavioral models (6b)', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

hold on;
plot(Bpr(2,Subjects(:,2)==1), Bpr(1,Subjects(:,2)==1), '.r', 'MarkerSize', 10);
plot(Bpr(2,Subjects(:,2)==2), Bpr(1,Subjects(:,2)==2), '.b', 'MarkerSize', 10);
xlim([min(Bpr(2,:))-1/10*range(Bpr(2,:)), max(Bpr(2,:))+1/10*range(Bpr(2,:))]);
ylim([min(Bpr(1,:))-1/10*range(Bpr(1,:)), max(Bpr(1,:))+1/10*range(Bpr(1,:))]);
axis square;
set(gca,'Box','On');
legend(grl, 'Location', 'NorthEast');
xlabel('\beta_1', 'FontSize', 12);
ylabel('\beta_2', 'FontSize', 12);
title('linear regression: Xp -> Zr', 'FontSize', 16);