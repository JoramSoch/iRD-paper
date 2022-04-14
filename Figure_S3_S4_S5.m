% NBD: Figures S3, S4, S5

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
DA_meas  =  'BA';

% search, if file does not exist
if ~exist(strcat('Figure_S5_S6_S7_',DA_meas,'.mat'), 'file')

    % preallocate results
    DAPs_beh = zeros(num_subj, 6, numel(ana_ids));
    DAPs_des = zeros(num_subj, 2, numel(ana_ids));
    i_subj   = [];

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
                    if strcmp(DA_meas,'DA')
                        DAPs_beh(i,:,k) = [PBM.DA, CRD.DA, NBD.DA];
                    elseif strcmp(DA_meas,'BA')
                        DAPs_beh(i,:,k) = [PBM.BA, CRD.BA, NBD.BA];
                    end;
                end;
                if strncmp(ana_ids{k},'3',1) || strncmp(ana_ids{k},'4',1)
                    DAPs_beh(i,:,k) = [PBM.r, CRD.r, NBD.r];
                end;

                % store results (design)
                if strncmp(ana_ids{k},'1',1) || strncmp(ana_ids{k},'3',1)
                    if strcmp(DA_meas,'DA')
                        DAPs_des(i,:,k) = NBD.DA_X;
                    elseif strcmp(DA_meas,'BA')
                        DAPs_des(i,:,k) = NBD.BA_X;
                    end;
                end;
                if strncmp(ana_ids{k},'2',1) || strncmp(ana_ids{k},'4',1)
                    DAPs_des(i,:,k) = NBD.r_X;
                end;

            end;

        end;

    end;

    % restrict results
    subj_ids = subj_ids(i_subj);
    subj_grp = Subjects(i_subj,2);
    DAPs_beh = DAPs_beh(i_subj,:,:);
    DAPs_des = DAPs_des(i_subj,:,:);

    % save results
    save(strcat('Figure_S5_S6_S7_',DA_meas,'.mat'), 'subj_ids', 'subj_grp', 'DAPs_beh', 'DAPs_des');
    
% load, if file does exist
else
    
    % load results
    load(strcat('Figure_S5_S6_S7_',DA_meas,'.mat'));
    
end;

% prepare plotting
X_labs = {'gain', 'loss'};
Z_labs = {'sbRD', 'dRD', 'iRD', '+gain', '+loss', '+both'};
Z_labs_LaTeX = {'$X \!\!\! \rightarrow \!\! Z$', '$Y_2 \!\! \rightarrow \!\! Z \quad$', '$Y_1 \!\! \rightarrow \!\! X \!\!\! \rightarrow \!\! Z$'};
y_sig  = [0.95, 0.95, 0.9, 0.9];

% plot results (PBM/CRD/NBD)
for i = 1:2
    % open figure
    figure('Name', sprintf('PBM/CRD/NBD (%s)', grp_labs{i}), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
    % show analyses
    for k = 1:numel(ana_ids)
        subplot(1, numel(ana_ids), k); hold on;
        DAPs = DAPs_beh(subj_grp==i,1:3,k);
        h = boxplot(DAPs, 'colors', 'rgb', 'symbol', '+k');
        set(h, 'LineWidth', 1);
        set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabelRotation', 0);
        set(gca, 'XTick', [1:3], 'XTickLabel', Z_labs_LaTeX, 'FontSize', 14);
        xlabel('decoding algorithm', 'FontSize', 16);
        if strncmp(ana_ids{k},'1',1) || strncmp(ana_ids{k},'2',1)
            plot([(1-0.5), (3+0.5)], [0.5, 0.5], '-k', 'LineWidth', 1);
            axis([(1-0.5), (3+0.5), 0.45, 1.05]);
            if strcmp(DA_meas,'DA')
                ylabel('decoding accuracy', 'FontSize', 16);
            elseif strcmp(DA_meas,'BA')
                ylabel('balanced accuracy', 'FontSize', 16);
            end;
        end;
        if strncmp(ana_ids{k},'3',1) || strncmp(ana_ids{k},'4',1)
            plot([(1-0.5), (3+0.5)], [0, 0], '-k', 'LineWidth', 1);
            axis([(1-0.5), (3+0.5), -0.1, +1.1]);
            ylabel('correlation coefficient', 'FontSize', 16);
        end;
        [h, p] = ttest(DAPs(:,2), DAPs(:,3));
        plot([(2+0.1), (3-0.1)], [y_sig(k), y_sig(k)], '-k', 'LineWidth', 1);
        if p < 0.001
            text(2+0.5, y_sig(k), '***', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Baseline');
        elseif p < 0.01
            text(2+0.5, y_sig(k), '**',  'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Baseline');
        elseif p < 0.05
            text(2+0.5, y_sig(k), '*',   'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Baseline');
        else
            text(2+0.5, y_sig(k), 'n.s.','FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
        end;
        title(sprintf('Analysis %s', ana_ids{k}(1)), 'FontSize', 20);
    end;
end;

% plot results (NBD variants)
for i = 1:2
    % open figure
    figure('Name', sprintf('NBD variants (%s)', grp_labs{i}), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
    % show analyses
    for k = 1:numel(ana_ids)
        subplot(1, numel(ana_ids), k); hold on;
        DAPs = DAPs_beh(subj_grp==i,3:6,k);
        h = boxplot(DAPs, 'colors', 'bmmr', 'symbol', '+k');
        set(h, 'LineWidth', 1);
      % set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabelRotation', 0);
        set(gca, 'XTick', [1:4], 'XTickLabel', Z_labs(3:6), 'FontSize', 14);
        xlabel('decoding algorithm', 'FontSize', 16);
        if strncmp(ana_ids{k},'1',1) || strncmp(ana_ids{k},'2',1)
            plot([(1-0.5), (4+0.5)], [0.5, 0.5], '-k', 'LineWidth', 1);
            axis([(1-0.5), (4+0.5), 0.45, 1.05]);
            if strcmp(DA_meas,'DA')
                ylabel('decoding accuracy', 'FontSize', 16);
            elseif strcmp(DA_meas,'BA')
                ylabel('balanced accuracy', 'FontSize', 16);
            end;
        end;
        if strncmp(ana_ids{k},'3',1) || strncmp(ana_ids{k},'4',1)
            plot([(1-0.5), (4+0.5)], [0, 0], '-k', 'LineWidth', 1);
            axis([(1-0.5), (4+0.5), -0.1, +1.1]);
            ylabel('correlation coefficient', 'FontSize', 16);
        end;
        [h, p] = ttest(DAPs(:,2), DAPs(:,3));
        plot([(2+0.1), (3-0.1)], [y_sig(k), y_sig(k)], '-k', 'LineWidth', 1);
        if p < 0.001
            text(2+0.5, y_sig(k), '***', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Baseline');
        elseif p < 0.01
            text(2+0.5, y_sig(k), '**',  'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Baseline');
        elseif p < 0.05
            text(2+0.5, y_sig(k), '*',   'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Baseline');
        else
            text(2+0.5, y_sig(k), 'n.s.','FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
        end;
        title(sprintf('Analysis %s', ana_ids{k}(1)), 'FontSize', 20);
    end;
end;

% plot results (NBD design)
y_sig = [0.95, 0.9, 0.95, 0.9];
for i = 1:2
    % open figure
    figure('Name', sprintf('NBD design (%s)', grp_labs{i}), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
    % show analyses
    for k = 1:numel(ana_ids)
        subplot(1, numel(ana_ids), k); hold on;
        DAPs = DAPs_des(subj_grp==i,:,k);
        h = boxplot(DAPs, 'colors', 'bb', 'symbol', '+k');
        set(h, 'LineWidth', 1);
      % set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabelRotation', 0);
        set(gca, 'XTick', [1:2], 'XTickLabel', X_labs, 'FontSize', 14);
        xlabel('decoded variable', 'FontSize', 16);
        if strncmp(ana_ids{k},'1',1) || strncmp(ana_ids{k},'3',1)
            plot([(1-0.5), (2+0.5)], [0.5, 0.5], '-k', 'LineWidth', 1);
            axis([(1-0.5), (2+0.5), 0.45, 1.05]);
            if strcmp(DA_meas,'DA')
                ylabel('decoding accuracy', 'FontSize', 16);
            elseif strcmp(DA_meas,'BA')
                ylabel('balanced accuracy', 'FontSize', 16);
            end;
        end;
        if strncmp(ana_ids{k},'2',1) || strncmp(ana_ids{k},'4',1)
            plot([(1-0.5), (2+0.5)], [0, 0], '-k', 'LineWidth', 1);
            axis([(1-0.5), (2+0.5), -0.1, +1.1]);
            ylabel('correlation coefficient', 'FontSize', 16);
        end;
        [h, p] = ttest(DAPs(:,1), DAPs(:,2));
        plot([(1+0.1), (2-0.1)], [y_sig(k), y_sig(k)], '-k', 'LineWidth', 1);
        if p < 0.001
            text(1+0.5, y_sig(k), '***', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Baseline');
        elseif p < 0.01
            text(1+0.5, y_sig(k), '**',  'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Baseline');
        elseif p < 0.05
            text(1+0.5, y_sig(k), '*',   'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Baseline');
        else
            text(1+0.5, y_sig(k), 'n.s.','FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
        end;
        title(sprintf('Analysis %s', ana_ids{k}(1)), 'FontSize', 20);
    end;
end;