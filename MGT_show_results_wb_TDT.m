function MGT_show_results_wb_TDT(data_dir, subj_id, MS_name, GLM_names, ana_id, plots)
% _
% Compare CRD vs. NBD decoding results (TDT)
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 03/03/2021, 12:46 (1/2/3/4d)


% things to do
thx2do   = [0 1 2];             % 0 = load, 1 = save, 2 = plot
subj_dir = strcat(data_dir,subj_id,'/glms/','glms-',MS_name);

% specify plots
if isempty(plots) || nargin < 6
    plots = 1;
end;
if plots == 0
    thx2do = [0 1];
end;


%%% Step 0: load everything %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(0,thx2do)

% load design information
load(strcat(subj_dir,'/','NBD_','XZ.mat'));

% load decoding results (PBM)
load(strcat(subj_dir,'/','NBD_','PBM-',ana_id,'.mat'));

% load decoding results (CRD)
if strcmp(ana_id,'1d') || strcmp(ana_id,'2d'), type_Z = 'SVC'; end;
if strcmp(ana_id,'3d') || strcmp(ana_id,'4d'), type_Z = 'SVR'; end;
TDT_Z = load(strcat(subj_dir,'/','glm-',GLM_names{1},'/','TDT_',type_Z,'_fav','/','results.mat')); % close
TDT_Z = TDT_Z.results;

% load decoding results (NBD)
if strcmp(ana_id,'1d') || strcmp(ana_id,'3d'), type_X = 'SVC'; end;
if strcmp(ana_id,'2d') || strcmp(ana_id,'4d'), type_X = 'SVR'; end;
TDT_X1 = load(strcat(subj_dir,'/','glm-',GLM_names{2},'/','TDT_',type_X,'_gain','/','results.mat')); % close
TDT_X2 = load(strcat(subj_dir,'/','glm-',GLM_names{2},'/','TDT_',type_X,'_loss','/','results.mat')); % close
TDT_X1 = TDT_X1.results;
TDT_X2 = TDT_X2.results;

end;


%%% Step 1: calculate accuracies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(1,thx2do)

%%% Intro: behavioral responses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analysis 1/3: categorical conditions
if strncmp(ana_id,'1',1) || strncmp(ana_id,'3',1)
    X = Xc;
end;
% Analysis 2/4: parametric regressors
if strncmp(ana_id,'2',1) || strncmp(ana_id,'4',1)
    X = Xp;
end;
% Analysis 1/2: categorical responses
if strncmp(ana_id,'1',1) || strncmp(ana_id,'2',1)
    Z = Zc;
end;
% Analysis 3/4: parametric responses
if strncmp(ana_id,'3',1) || strncmp(ana_id,'4',1)
    Z = Zp;
end;


%%% PBM: decode Z from X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear PBM

% Analysis 1: categorical conditions -> categorical responses
if strncmp(ana_id,'1',1)
    PBM.Z = Zcc;
end;
% Analysis 2: parametric regressors -> categorical responses
if strncmp(ana_id,'2',1)
    PBM.Z = Zpc;
end;
%  Analysis 1/2: categorical responses
if strncmp(ana_id,'1',1)
   [PBM.DA, PBM.BA] = NBD_calc_DA(Z(:,2), PBM.Z(:,2), iZ);
end;
if strncmp(ana_id,'2',1)
   [PBM.DA, PBM.BA] = NBD_calc_DA(Z(:,2), PBM.Z, iZ);
end;
if strncmp(ana_id,'1',1) || strncmp(ana_id,'2',1)
   [p, PBM.DA_CI] = binofit(PBM.DA*sum(iZ), sum(iZ));
   [p, PBM.BA_CI] = binofit(floor(PBM.BA*sum(iZ)), sum(iZ));
    PBM.DA_CI = PBM.DA_CI';
    PBM.BA_CI = PBM.BA_CI';
end;
% Analysis 3: categorical conditions -> parametric responses
if strncmp(ana_id,'3',1)
    PBM.Z = Zcp;
end;
% Analysis 4: parametric regressors -> parametric responses
if strncmp(ana_id,'4',1)
    PBM.Z = Zpp;
end;
% Analysis 3/4: parametric responses
if strncmp(ana_id,'3',1) || strncmp(ana_id,'4',1)
   [PBM.r, PBM.p, PBM.MSE, PBM.MAE, PBM.slo, PBM.int] = NBD_calc_DP(Z, PBM.Z, iZ);
   [R, P, RL, RU] = corrcoef(Z(iZ), PBM.Z(iZ));
    PBM.r_lu = [RL(1,2); RU(1,2)];
end;


%%% CRD: decode Z from Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear CRD

% Analysis 1/2: measured data -> categorical responses
if strncmp(ana_id,'1',1) || strncmp(ana_id,'2',1)
    CRD.Z = zeros(size(Zc,1),1);
    CRD.Z(iZ) = TDT_Z.predicted_labels.output{1}-1;
   [CRD.DA, CRD.BA] = NBD_calc_DA(Z(:,2), CRD.Z, iZ);
   [p, CRD.DA_CI] = binofit(round(CRD.DA*sum(iZ)), sum(iZ));
   [p, CRD.BA_CI] = binofit(floor(CRD.BA*sum(iZ)), sum(iZ));
    CRD.DA_CI = CRD.DA_CI';
    CRD.BA_CI = CRD.BA_CI';
end;
% Analysis 3/4: measured data -> parametric responses
if strncmp(ana_id,'3',1) || strncmp(ana_id,'4',1)
    CRD.Z = zeros(size(Zp,1),1);
    CRD.Z(iZ) = TDT_Z.predicted_labels.output{1};
   [CRD.r, CRD.p, CRD.MSE, CRD.MAE, CRD.slo, CRD.int] = NBD_calc_DP(Z, CRD.Z, iZ);
   [R, P, RL, RU] = corrcoef(Z(iZ), CRD.Z(iZ));
    CRD.r_lu = [RL(1,2); RU(1,2)];
end;


%%% NBD: decode Z from Y via X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear NBD

% Analysis 1/3: measured data -> categorical conditions
if strncmp(ana_id,'1',1) || strncmp(ana_id,'3',1)
    NBD.x = zeros(size(Xc,1),log2(size(Xc,2)));
    for k = 1:size(NBD.x,2)
        if k == 1, NBD.x(:,k) = 2-TDT_X1.predicted_labels.output{1}; end; % gain
        if k == 2, NBD.x(:,k) = 2-TDT_X2.predicted_labels.output{1}; end; % loss
        if k == 1, Xk = [sum(X(:,[1 2]),2), sum(X(:,[3 4]),2)]; end; % gain
        if k == 2, Xk = [sum(X(:,[1 3]),2), sum(X(:,[2 4]),2)]; end; % loss
       [NBD.DA_X(k), NBD.BA_X(k)] = NBD_calc_DA(Xk(:,1), NBD.x(:,k));
       [p, NBD.DA_CI_X(k,:)] = binofit(round(NBD.DA_X(k)*numel(Xk(:,1))), numel(Xk(:,1)));
       [p, NBD.BA_CI_X(k,:)] = binofit(floor(NBD.BA_X(k)*numel(Xk(:,1))), numel(Xk(:,1)));
    end;
    NBD.DA_CI_X = NBD.DA_CI_X';
    NBD.BA_CI_X = NBD.BA_CI_X';
    xA = sum(X(:,[1 2]),2);     % low gain
    xB = sum(X(:,[1 3]),2);     % low loss
    NBD.X = zeros(size(Xc,1),size(Xc,2),2^log2(size(Xc,2)));
    NBD.X(:,:,1) = NBD_marg2joint([NBD.x(:,1), NBD.x(:,2)]);
    NBD.X(:,:,2) = NBD_marg2joint([NBD.x(:,1), xB        ]);
    NBD.X(:,:,3) = NBD_marg2joint([xA,         NBD.x(:,2)]);
    NBD.X(:,:,4) = NBD_marg2joint([xA,         xB        ]);
end;
% Analysis 2/4: measured data -> parametric regressors
if strncmp(ana_id,'2',1) || strncmp(ana_id,'4',1)
    NBD.x = zeros(size(Xp,1),size(Xp,2)-1);
    for k = 1:size(NBD.x,2)
        if k == 1, NBD.x(:,k) = TDT_X1.predicted_labels.output{1}; end; % gain
        if k == 2, NBD.x(:,k) = TDT_X2.predicted_labels.output{1}; end; % loss
       [NBD.r_X(k), NBD.p_X(k), NBD.MSE_X(k), NBD.MAE_X(k), NBD.slo_X(k), NBD.int_X(k)] = NBD_calc_DP(X(:,k), NBD.x(:,k));
       [R, P, RL, RU]   = corrcoef(X(:,k), NBD.x(:,k));
        NBD.r_lu_X(:,k) = [RL(1,2); RU(1,2)];
    end;
    NBD.X = zeros(size(Xp,1),size(Xp,2),2^(size(Xp,2)-1));
    NBD.X(:,:,1) = [NBD.x(:,1), NBD.x(:,2), X(:,3)];
    NBD.X(:,:,2) = [NBD.x(:,1), X(:,2),     X(:,3)];
    NBD.X(:,:,3) = [X(:,1),     NBD.x(:,2), X(:,3)];
    NBD.X(:,:,4) = [X(:,1),     X(:,2),     X(:,3)];
end;
% Analysis 1: categorical conditions -> categorical responses
if strncmp(ana_id,'1',1)
    NBD.Z = zeros(size(Zc,1),size(Zc,2),size(NBD.X,3));
    for f = 1:size(NBD.X,3)
        for g = 1:S
            NBD.Z(i2{g},:,f) = NBD.X(i2{g},:,f) * Bcc(:,:,g);
        end;
       [NBD.DA(f), NBD.BA(f)] = NBD_calc_DA(Z(:,2), NBD.Z(:,2,f), iZ);
       [p, NBD.DA_CI(f,:)] = binofit(round(NBD.DA(f)*sum(iZ)), sum(iZ));
       [p, NBD.BA_CI(f,:)] = binofit(floor(NBD.BA(f)*sum(iZ)), sum(iZ));
    end;
    NBD.DA_CI = NBD.DA_CI';
    NBD.BA_CI = NBD.BA_CI';
end;
% Analysis 2: parametric regressors -> categorical responses
if strncmp(ana_id,'2',1)
    NBD.Z = zeros(size(Zc,1),size(NBD.X,3));
    for f = 1:size(NBD.X,3)
        for g = 1:S
            NBD.Z(i2{g},f) = 1./(1 + exp(-NBD.X(i2{g},:,f)*Bpc(:,g)));
        end;
       [NBD.DA(f), NBD.BA(f)] = NBD_calc_DA(Z(:,2), NBD.Z(:,f), iZ);
       [p, NBD.DA_CI(f,:)] = binofit(round(NBD.DA(f)*sum(iZ)), sum(iZ));
       [p, NBD.BA_CI(f,:)] = binofit(floor(NBD.BA(f)*sum(iZ)), sum(iZ));
    end;
    NBD.DA_CI = NBD.DA_CI';
    NBD.BA_CI = NBD.BA_CI';
end;
% Analysis 3: categorical conditions -> parametric responses
if strncmp(ana_id,'3',1)
    NBD.Z = zeros(size(Zp,1),size(NBD.X,3));
    for f = 1:size(NBD.X,3)
        for g = 1:S
            NBD.Z(i2{g},f) = NBD.X(i2{g},:,f) * Bcp(:,g);
        end;
       [NBD.r(f), NBD.p(f), NBD.MSE(f), NBD.MAE(f), NBD.slo(f), NBD.int(f)] = NBD_calc_DP(Z, NBD.Z(:,f), iZ);
       [R, P, RL, RU] = corrcoef(Z(iZ), NBD.Z(iZ,f));
        NBD.r_lu(:,f) = [RL(1,2); RU(1,2)];
    end;
end;
% Analysis 4: parametric regressors -> parametric responses
if strncmp(ana_id,'4',1)
    NBD.Z = zeros(size(Zp,1),size(NBD.X,3));
    for f = 1:size(NBD.X,3)
        for g = 1:S
            NBD.Z(i2{g},f) = NBD.X(i2{g},:,f) * Bpp(:,g);
        end;
       [NBD.r(f), NBD.p(f), NBD.MSE(f), NBD.MAE(f), NBD.slo(f), NBD.int(f)] = NBD_calc_DP(Z, NBD.Z(:,f), iZ);
       [R, P, RL, RU] = corrcoef(Z(iZ), NBD.Z(iZ,f));
        NBD.r_lu(:,f) = [RL(1,2); RU(1,2)];
    end;
end;

% save results
clear p R P RL RU
filename = strcat(subj_dir,'/','Res_NBD-',ana_id,'.mat');
save(filename, 'X', 'Z', 'iZ', 'PBM', 'CRD', 'NBD');

end;


%%% Step 2: display accuracies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(2,thx2do)

% close all

% load results
filename = strcat(subj_dir,'/','Res_NBD-',ana_id,'.mat');
load(filename);

% plot measured data
if plots > 1

figure('Name', 'NBD: X & Z', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

% behavior
subplot(1,5,1);
imagesc(Z);
if strncmp(ana_id,'1',1) || strncmp(ana_id,'2',1), caxis([0 1]); end;
if strncmp(ana_id,'3',1) || strncmp(ana_id,'4',1), caxis([0 4]); end;
colorbar;
set(gca,'XTick',[1:size(Z,2)]);
xlabel('response', 'FontSize', 12);
ylabel('trial', 'FontSize', 12);
title('behavioral responses Z', 'FontSize', 16);

% design
subplot(1,5,[2:4]);
imagesc(X);
if strncmp(ana_id,'1',1) || strncmp(ana_id,'3',1), caxis([0 1]); end;
if strncmp(ana_id,'2',1) || strncmp(ana_id,'4',1), caxis([min(min(X)), max(max(X))]); end;
colorbar;
set(gca,'XTick',[1:size(X,2)]);
xlabel('condition', 'FontSize', 12);
ylabel('trial', 'FontSize', 12);
title('experimental design X', 'FontSize', 16);

% parameters
subplot(1,5,5);
if strncmp(ana_id,'1',1), B = mean(Bcc,3); end;
if strncmp(ana_id,'2',1), B = mean(Bpc,2); end;
if strncmp(ana_id,'3',1), B = mean(Bcp,2); end;
if strncmp(ana_id,'4',1), B = mean(Bpp,2); end;
imagesc(B);
if strncmp(ana_id,'1',1), caxis([0 1]); end;
colorbar;
set(gca,'XTick',[1:size(B,2)]);
set(gca,'YTick',[1:size(B,1)]);
xlabel('response', 'FontSize', 12);
ylabel('condition', 'FontSize', 12);
title('mapping parameters P', 'FontSize', 16);

end;

% plot decoding accuracies/precisions
if plots > 0

h = figure('Name', 'NBD: DAs/DPs', 'Color', [1 1 1], 'Position', [1280 -120 1600 900]);
Z_labs = {'PBM', 'CRD', 'NBD', 'NBD+gain', 'NBD+loss', 'NBD+both'};
X_labs = {'gain', 'loss'};

% behavior
if strncmp(ana_id,'1',1) || strncmp(ana_id,'2',1)
    subplot(1,3,[1:2]); hold on;
    DAs = [PBM.BA, CRD.BA, NBD.BA];
    CIs = [PBM.BA_CI, CRD.BA_CI, NBD.BA_CI];
    bar([1:numel(DAs)], DAs, 'b');
    errorbar([1:numel(DAs)], DAs, DAs-CIs(1,:), CIs(2,:)-DAs, '.k', 'LineWidth', 2, 'CapSize', 15);
    axis([(1-1), (numel(DAs)+1), -0.01, +1.01]);
    set(gca,'Box','On');
    set(gca, 'XTick', [1:numel(DAs)], 'XTickLabel', Z_labs);
    xlabel('decoding algorithm', 'FontSize', 12);
    ylabel('balanced accuracy', 'FontSize', 12);
    title('decoding accuracies (behavior)', 'FontSize', 16);
end;
if strncmp(ana_id,'3',1) || strncmp(ana_id,'4',1)
    subplot(1,3,[1:2]); hold on;
    DPs = [PBM.r, CRD.r, NBD.r];
    CIs = [PBM.r_lu, CRD.r_lu, NBD.r_lu];
    bar([1:numel(DPs)], DPs, 'b');
    errorbar([1:numel(DPs)], DPs, DPs-CIs(1,:), CIs(2,:)-DPs, '.k', 'LineWidth', 2, 'CapSize', 15);
    axis([(1-1), (numel(DPs)+1), -0.01, +1.01]);
    set(gca,'Box','On');
    set(gca, 'XTick', [1:numel(DPs)], 'XTickLabel', Z_labs);
    xlabel('decoding algorithm', 'FontSize', 12);
    ylabel('correlation coefficient', 'FontSize', 12);
    title('decoding precisions (behavior)', 'FontSize', 16);
end;
% design
if strncmp(ana_id,'1',1) || strncmp(ana_id,'3',1)
    subplot(1,3,3); hold on;
    DAs = NBD.BA_X;
    CIs = NBD.BA_CI_X;
    bar([1:numel(DAs)], DAs, 'r');
    axis([(1-1), (numel(DAs)+1), -0.01, +1.01]);
    errorbar([1:numel(DAs)], DAs, DAs-CIs(1,:), CIs(2,:)-DAs, '.k', 'LineWidth', 2, 'CapSize', 15);
    set(gca,'Box','On');
    set(gca, 'XTick', [1:numel(DAs)], 'XTickLabel', X_labs);
    xlabel('experimental variable', 'FontSize', 12);
    ylabel('balanced accuracy', 'FontSize', 12);
    title('decoding accuracies (design)', 'FontSize', 16);
end;
if strncmp(ana_id,'2',1) || strncmp(ana_id,'4',1)
    subplot(1,3,3); hold on;
    DPs = NBD.r_X;
    CIs = NBD.r_lu_X;
    bar([1:numel(DPs)], DPs, 'r');
    axis([(1-1), (numel(DPs)+1), -0.01, +1.01]);
    errorbar([1:numel(DPs)], DPs, DPs-CIs(1,:), CIs(2,:)-DPs, '.k', 'LineWidth', 2, 'CapSize', 15);
    set(gca,'Box','On');
    set(gca, 'XTick', [1:numel(DPs)], 'XTickLabel', X_labs);
    xlabel('experimental variable', 'FontSize', 12);
    ylabel('correlation coefficient', 'FontSize', 12);
    title('decoding precisions (design)', 'FontSize', 16);
end;

filename = strcat(subj_dir,'/','Res_NBD-',ana_id,'.png');
saveas(h, filename);

end;

end;