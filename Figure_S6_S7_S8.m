% NBD: Figures S6, S7, S8

clear
close all

% load subjects
tool_dir = strcat(pwd,'/');
load(strcat(tool_dir,'extract_subjects.mat'));
load(strcat(tool_dir,'project_directories.mat'));
load(strcat(tool_dir,'MS_files/MS_MGT1.mat'));

% specify analyses
num_subj = numel(subj_ids);
ana_ids  = {'1d', '2d', '3d', '4d'};
grp_labs = {'equal range', 'equal indifference'};
DA_meas  =  'DA';

% specify statistics
chlvl = 0.5;                    % chance level of first-level decoding
alpha = 0.05;                   % significance level for first-level performance
beta  = 1;                      % statistical power for first-level performance
pdens = 0.9;                    % posterior density level for second-level inference

% search, if file does not exist
if ~exist(strcat('Figure_S8_S9_S10_',DA_meas,'.mat'), 'file')

    % preallocate results
    ps_beh = zeros(num_subj, 6, numel(ana_ids));
    ps_des = zeros(num_subj, 2, numel(ana_ids));
    i_subj = [];

    % for all subjects
    for i = 1:num_subj

        % specify subject directory
        subj_dir = strcat(data_dir,subj_ids{i},'/glms/','glms-',MS_name);
        if exist(subj_dir,'dir')
            i_subj = [i_subj, i];

            % for all analyses
            for k = 1:numel(ana_ids)

                % specify results file
                filename = strcat(subj_dir,'/','Res_NBD-',ana_ids{k},'.mat');
                load(filename);

                % store results (behavior)
                if strncmp(ana_ids{k},'1',1) || strncmp(ana_ids{k},'2',1)
                    % Explanation: p-value is set to 0, if CI does not
                    % contain chance level, and to 1 otherwise.
                    if strcmp(DA_meas,'DA')
                        ps_beh(i,:,k) = double([PBM.DA_CI(1,:), CRD.DA_CI(1,:), NBD.DA_CI(1,:)] < chlvl);
                    elseif strcmp(DA_meas,'BA')
                        ps_beh(i,:,k) = double([PBM.BA_CI(1,:), CRD.BA_CI(1,:), NBD.BA_CI(1,:)] < chlvl);
                    end;
                end;
                if strncmp(ana_ids{k},'3',1) || strncmp(ana_ids{k},'4',1)
                    ps_beh(i,:,k) = [PBM.p, CRD.p, NBD.p];
                end;

                % store results (design)
                if strncmp(ana_ids{k},'1',1) || strncmp(ana_ids{k},'3',1)
                    % Explanation: p-value is set to 0, if CI does not
                    % contain chance level, and to 1 otherwise.
                    if strcmp(DA_meas,'DA')
                        ps_des(i,:,k) = double(NBD.DA_CI_X(1,:) < chlvl);
                    elseif strcmp(DA_meas,'BA')
                        ps_des(i,:,k) = double(NBD.BA_CI_X(1,:) < chlvl);
                    end;
                end;
                if strncmp(ana_ids{k},'2',1) || strncmp(ana_ids{k},'4',1)
                    ps_des(i,:,k) = NBD.p_X;
                end;

            end;

        end;

    end;

    % statistically test
    hs_beh = double(ps_beh < alpha);
    hs_des = double(ps_des < alpha);

    % restrict results
    subj_ids = subj_ids(i_subj);
    subj_grp = Subjects(i_subj,2);
    hs_beh   = hs_beh(i_subj,:,:);
    hs_des   = hs_des(i_subj,:,:);
    
    % save results
    save(strcat('Figure_S8_S9_S10_',DA_meas,'.mat'), 'subj_ids', 'subj_grp', 'hs_beh', 'hs_des');
    
% load, if file does exist
else
    
    % load results
    load(strcat('Figure_S8_S9_S10_',DA_meas,'.mat'));
    
end;

% prepare plotting
n      = [sum(subj_grp==1), sum(subj_grp==2), numel(subj_ids)];
X_labs = {'gain', 'loss'};
Z_labs = {'sbRD', 'dRD', 'iRD', '+gain', '+loss', '+both'};
Z_labs_LaTeX = {'$X \!\!\! \rightarrow \!\! Z$', '$Y_2 \!\! \rightarrow \!\! Z \quad$', '$Y_1 \!\! \rightarrow \!\! X \!\!\! \rightarrow \!\! Z$'};

% plot results (PBM/CRD/NBD)
for i = 1:2
    % open figure
    figure('Name', sprintf('PBM/CRD/NBD (%s)', grp_labs{i}), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
    cols = 'rgb';
    % show analyses
    for k = 1:numel(ana_ids)
        subplot(1, numel(ana_ids), k); hold on;
        k_i    = sum(hs_beh(subj_grp==i,1:3,k));
        g_MAP  = zeros(1,numel(k_i));
        g_hpdi = zeros(2,numel(k_i));
        for l = 1:numel(k_i)
            g_MAP(l)    = bayesprev_map(k_i(l), n(i), alpha, beta);
            g_hpdi(:,l) = bayesprev_hpdi(pdens, k_i(l), n(i), alpha, beta)';
            bar(l, g_MAP(l), cols(l));
            errorbar(l, g_MAP(l), (g_MAP(l)-g_hpdi(1,l)), (g_hpdi(2,l)-g_MAP(l)), '.k', 'LineWidth', 2, 'CapSize', 20);
        end;
        axis([(1-0.6), (numel(k_i)+0.6), -0.05, +1.05]);
        set(gca, 'Box', 'On');
        set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabelRotation', 0);
        set(gca, 'XTick', [1:3], 'XTickLabel', Z_labs_LaTeX, 'FontSize', 14);
        xlabel('decoding algorithm', 'FontSize', 16);
        ylabel('effect prevalence', 'FontSize', 16);
      % ylabel('effect prevalence [MAP & hpdi]', 'FontSize', 16);
        title(sprintf('Analysis %s', ana_ids{k}(1)), 'FontSize', 20);
    end;
end;

% plot results (NBD variants)
for i = 1:2
    % open figure
    figure('Name', sprintf('NBD variants (%s)', grp_labs{i}), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
    cols = 'bmmr';
    % show analyses
    for k = 1:numel(ana_ids)
        subplot(1, numel(ana_ids), k); hold on;
        k_i    = sum(hs_beh(subj_grp==i,3:6,k));
        g_MAP  = zeros(1,numel(k_i));
        g_hpdi = zeros(2,numel(k_i));
        for l = 1:numel(k_i)
            g_MAP(l)    = bayesprev_map(k_i(l), n(i), alpha, beta);
            g_hpdi(:,l) = bayesprev_hpdi(pdens, k_i(l), n(i), alpha, beta)';
            bar(l, g_MAP(l), cols(l));
            errorbar(l, g_MAP(l), (g_MAP(l)-g_hpdi(1,l)), (g_hpdi(2,l)-g_MAP(l)), '.k', 'LineWidth', 2, 'CapSize', 20);
        end;
        axis([(1-0.6), (numel(k_i)+0.6), -0.05, +1.05]);
        set(gca, 'Box', 'On');
        set(gca, 'XTick', [1:4], 'XTickLabel', Z_labs(3:6), 'FontSize', 14);
        xlabel('decoding algorithm', 'FontSize', 16);
        ylabel('effect prevalence', 'FontSize', 16);
      % ylabel('effect prevalence [MAP & hpdi]', 'FontSize', 16);
        title(sprintf('Analysis %s', ana_ids{k}(1)), 'FontSize', 20);
    end;
end;

% plot results (NBD design)
for i = 1:2
    % open figure
    figure('Name', sprintf('NBD design (%s)', grp_labs{i}), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
    cols = 'bb';
    % show analyses
    for k = 1:numel(ana_ids)
        subplot(1, numel(ana_ids), k); hold on;
        k_i    = sum(hs_des(subj_grp==i,:,k));
        g_MAP  = zeros(1,numel(k_i));
        g_hpdi = zeros(2,numel(k_i));
        for l = 1:numel(k_i)
            g_MAP(l)    = bayesprev_map(k_i(l), n(i), alpha, beta);
            g_hpdi(:,l) = bayesprev_hpdi(pdens, k_i(l), n(i), alpha, beta)';
            bar(l, g_MAP(l), cols(l));
            errorbar(l, g_MAP(l), (g_MAP(l)-g_hpdi(1,l)), (g_hpdi(2,l)-g_MAP(l)), '.k', 'LineWidth', 2, 'CapSize', 20);
        end;
        axis([(1-1), (numel(k_i)+1), -0.05, +1.05]);
        set(gca,'Box','On');
        set(gca, 'XTick', [1:2], 'XTickLabel', X_labs);
        xlabel('decoded variable', 'FontSize', 16);
        ylabel('effect prevalence', 'FontSize', 16);
      % ylabel('effect prevalence [MAP & hpdi]', 'FontSize', 16);
        title(sprintf('Analysis %s', ana_ids{k}(1)), 'FontSize', 20);
    end;
end;