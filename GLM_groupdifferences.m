% 1. LOAD CONFIGURATION AND SPREADSHEET
p = opm_study_config_babyeinsteins_NNMF();

% Pointing to the new Results_1 folder where the updated spreadsheet lives
table_file = '/path/to/nnmf_results.xlsx'; 

% Load table
data = readtable(table_file);

% 2. PREPARE VARIABLES FOR STANDARDIZED BETAS
data.dx = categorical(strtrim(data.dx), {'TDC', 'ASD'});
data.sex = categorical(strtrim(data.sex)); 

% Ensure age is numeric and Z-SCORE IT
if ~isnumeric(data.age)
    data.age = str2double(data.age);
end
data.age = (data.age - mean(data.age, 'omitnan')) / std(data.age, 'omitnan');

fprintf('Running Doherty Robust GLM (Controlling for Age & Sex) on RAW ICA data...\n\n');

% 3. IDENTIFY ALL RAW ICA COLUMNS 
all_vars = data.Properties.VariableNames;

% We search specifically for the RAW columns for valid cross-sectional group comparisons
raw_cols = all_vars(contains(all_vars, '_Raw_70_ica_pca')); 

num_components = length(raw_cols);
if num_components == 0
    error('No columns containing "_Raw_70_ica_pca" were found in the table. Check column names!');
end

fprintf('Detected %d Total Raw Components.\n', num_components);

% 4. INITIALIZE STORAGE ARRAYS
std_beta_group = zeros(num_components, 1);
ci_lower       = zeros(num_components, 1);
ci_upper       = zeros(num_components, 1);
p_vals         = zeros(num_components, 1);
outliers_drop  = zeros(num_components, 1);
p_age          = zeros(num_components, 1);
p_sex          = zeros(num_components, 1);
models         = cell(num_components, 1); % Storage to save models for plotting later

% 5. RUN ROBUST LINEAR MODELS
for i = 1:num_components
    col_name = raw_cols{i};
    
    % Z-SCORE THE NETWORK STRENGTH (Required for standardized betas)
    y_raw = data.(col_name);
    data.y_std = (y_raw - mean(y_raw, 'omitnan')) / std(y_raw, 'omitnan');
    
    % Define the formula
    formula = 'y_std ~ dx + sex + age';
    
    % --- STEP 1: INITIAL FIT ---
    mdl_initial = fitlm(data, formula);
    
    % --- STEP 2: COOK's DISTANCE OUTLIER REJECTION ---
    cooksD = mdl_initial.Diagnostics.CooksDistance;
    mean_cooks = mean(cooksD, 'omitnan');
    
    % Doherty Rules: > 3x mean OR absolute value > 0.5
    outlier_idx = (cooksD > (3 * mean_cooks)) | (abs(cooksD) > 0.5);
    outliers_drop(i) = sum(outlier_idx);
    
    % Create a clean dataset without the outliers
    data_clean = data;
    data_clean(outlier_idx, :) = [];
    
    % --- STEP 3: FINAL ROBUST FIT ---
    mdl_final = fitlm(data_clean, formula);
    
    % Store the final model so we can plot it later
    models{i} = mdl_final;
    
    % Extract the Standardized Beta for Group (dx_ASD)
    std_beta_group(i) = mdl_final.Coefficients{'dx_ASD', 'Estimate'};
    p_vals(i)         = mdl_final.Coefficients{'dx_ASD', 'pValue'};
    
    % Extract 95% Confidence Intervals
    ci = coefCI(mdl_final);
    dx_idx = strcmp(mdl_final.CoefficientNames, 'dx_ASD');
    ci_lower(i) = ci(dx_idx, 1);
    ci_upper(i) = ci(dx_idx, 2);
    
    % Extract Covariate p-values
    sex_coef_name = mdl_final.CoefficientNames(contains(mdl_final.CoefficientNames, 'sex_'));
    if ~isempty(sex_coef_name)
        p_sex(i) = mdl_final.Coefficients{sex_coef_name{1}, 'pValue'};
    else
        p_sex(i) = NaN; 
    end
    p_age(i) = mdl_final.Coefficients{'age', 'pValue'};
end

% 6. SEPARATE FDR (BENJAMINI-HOCHBERG) CORRECTION 
alpha_thresh = 0.05;
p_fdr = zeros(num_components, 1);
tests_applied = zeros(num_components, 1);

act_idx = find(contains(raw_cols, 'Activity'));
conn_idx = find(contains(raw_cols, 'Connectivity'));

% FDR for Activity
if ~isempty(act_idx)
    n_act = length(act_idx);
    [p_sorted_act, sort_idx_act] = sort(p_vals(act_idx));
    padj_act = zeros(n_act, 1);
    padj_act(n_act) = p_sorted_act(n_act);
    for j = n_act-1:-1:1
        padj_act(j) = min(padj_act(j+1), p_sorted_act(j) * n_act / j);
    end
    fdr_mapped_act = zeros(n_act, 1);
    fdr_mapped_act(sort_idx_act) = padj_act;
    
    p_fdr(act_idx) = fdr_mapped_act;
    tests_applied(act_idx) = n_act;
end

% FDR for Connectivity
if ~isempty(conn_idx)
    n_conn = length(conn_idx);
    [p_sorted_conn, sort_idx_conn] = sort(p_vals(conn_idx));
    padj_conn = zeros(n_conn, 1);
    padj_conn(n_conn) = p_sorted_conn(n_conn);
    for j = n_conn-1:-1:1
        padj_conn(j) = min(padj_conn(j+1), p_sorted_conn(j) * n_conn / j);
    end
    fdr_mapped_conn = zeros(n_conn, 1);
    fdr_mapped_conn(sort_idx_conn) = padj_conn;
    
    p_fdr(conn_idx) = fdr_mapped_conn;
    tests_applied(conn_idx) = n_conn;
end

p_fdr = min(p_fdr, 1.0); % Cap at 1.0
fdr_sig = categorical(p_fdr <= alpha_thresh, [false, true], {'No', 'Yes*'});

% 7. DISPLAY FINAL SUMMARY TABLE
summary_table = table(raw_cols', std_beta_group, ci_lower, ci_upper, p_vals, p_fdr, fdr_sig, tests_applied, outliers_drop, p_age, p_sex, ...
    'VariableNames', {'Network', 'Std_Beta_ASD', 'CI_95_Low', 'CI_95_High', 'p_Group_Raw', 'p_FDR', 'FDR_Sig', 'Tests_Corrected_For', 'Outliers_Removed', 'p_Age', 'p_Sex'});

disp('========================================================================================');
disp('          DOHERTY ROBUST GLM RESULTS (Network ~ dx + sex + age)                         ');
disp('========================================================================================');
disp(summary_table);

% Save results
report_path = fullfile('/path/to/results/output/', 'Stats_Robust_GLM_GroupDiff_Raw_FDR.csv');
writetable(summary_table, report_path);
fprintf('\nSaved complete statistical report to: %s\n', report_path);

% 8. GENERATE PUBLICATION-READY VIOLIN PLOTS (SIGNIFICANT ONLY)
fprintf('\nGenerating Publication-Ready Violin Plots for SIGNIFICANT results...\n');
plot_dir = fullfile('/path/to/results/output/', 'GLM_Publication_Plots_Raw_FDR');

if ~exist(plot_dir, 'dir')
    mkdir(plot_dir); 
end

% Find only the models that survived FDR correction
sig_indices = find(p_fdr <= alpha_thresh);

if isempty(sig_indices)
    fprintf('No significant results found after FDR correction. No plots generated.\n');
else
    for idx = 1:length(sig_indices)
        i = sig_indices(idx);
        mdl = models{i};
        raw_name = raw_cols{i};
        
        % --- A. FORMAT A PUBLICATION-READY TITLE ---
        % Convert standard prefixes to readable frequency bands
        if contains(raw_name, 'x4_7')
            freq_label = 'Theta (4-7 Hz)';
        elseif contains(raw_name, 'x7_13')
            freq_label = 'Alpha (7-13 Hz)';
        elseif contains(raw_name, 'x13_30')
            freq_label = 'Beta (13-30 Hz)';
        elseif contains(raw_name, 'x30_55')
            freq_label = 'Low Gamma (30-55 Hz)';
        elseif contains(raw_name, 'x65_80')
            freq_label = 'High Gamma (65-80 Hz)';
        else
            freq_label = 'Network';
        end
        
        % Identify Modality
        if contains(raw_name, 'Connectivity')
            mod_label = 'Connectivity Network';
            tokens = regexp(raw_name, 'Connectivity_(\d+)', 'tokens');
        else
            mod_label = 'Activity Component';
            tokens = regexp(raw_name, 'Activity_(\d+)', 'tokens');
        end
        
        % Extract the network number
        if ~isempty(tokens)
            net_num = tokens{1}{1};
        else
            net_num = 'X';
        end
        
        clean_title = sprintf('%s %s %s', freq_label, mod_label, net_num);
        
        % --- B. EXTRACT PARTIAL RESIDUALS (Adjusted Scores) ---
        % We isolate the raw residuals and add back the effect of Group (dx)
        % This mathematically removes the variance caused by Age and Sex
        b0 = mdl.Coefficients{'(Intercept)', 'Estimate'};
        b_asd = mdl.Coefficients{'dx_ASD', 'Estimate'};
        is_asd = (mdl.Variables.dx == 'ASD');
        
        y_adj = mdl.Residuals.Raw + b0 + (b_asd .* is_asd);
        
        y_tdc = y_adj(mdl.Variables.dx == 'TDC');
        y_asd = y_adj(mdl.Variables.dx == 'ASD');
        
        % --- C. DRAW THE CUSTOM VIOLIN PLOT ---
        fig = figure('Name', clean_title, 'Color', 'w', 'Position', [200, 200, 500, 600], 'Visible', 'off');
        hold on;
        
        % Define professional colors
        color_TDC = [0 0.4470 0.7410]; % Blue
        color_ASD = [0.8500 0.3250 0.0980]; % Orange
        
        % 1. VIOLIN PATCHES (Using Kernel Density Estimation)
        % TDC Violin
        [f_tdc, xi_tdc] = ksdensity(y_tdc);
        f_tdc = 0.35 * (f_tdc / max(f_tdc)); % Scale width
        patch([1 - f_tdc, fliplr(1 + f_tdc)], [xi_tdc, fliplr(xi_tdc)], color_TDC, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % ASD Violin
        [f_asd, xi_asd] = ksdensity(y_asd);
        f_asd = 0.35 * (f_asd / max(f_asd)); % Scale width
        patch([2 - f_asd, fliplr(2 + f_asd)], [xi_asd, fliplr(xi_asd)], color_ASD, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % 2. BOXPLOTS (Clear, empty boxes showing quartiles)
        boxchart(ones(size(y_tdc)), y_tdc, 'BoxFaceColor', 'k', 'BoxFaceAlpha', 0, 'BoxWidth', 0.1, 'WhiskerLineColor', 'k', 'MarkerStyle', 'none', 'LineWidth', 1.2);
        boxchart(2*ones(size(y_asd)), y_asd, 'BoxFaceColor', 'k', 'BoxFaceAlpha', 0, 'BoxWidth', 0.1, 'WhiskerLineColor', 'k', 'MarkerStyle', 'none', 'LineWidth', 1.2);
        
        % 3. SWARMCHARTS (Individual Data Dots)
        swarmchart(ones(size(y_tdc)), y_tdc, 30, color_TDC, 'filled', 'MarkerFaceAlpha', 0.8, 'XJitterWidth', 0.15);
        swarmchart(2*ones(size(y_asd)), y_asd, 30, color_ASD, 'filled', 'MarkerFaceAlpha', 0.8, 'XJitterWidth', 0.15);
        
        % --- D. FORMATTING & AESTHETICS ---
        set(gca, 'XTick', [1, 2], 'XTickLabel', {'TDC', 'ASD'}, 'FontSize', 13, 'LineWidth', 1.5, 'Box', 'off');
        xlim([0.4, 2.6]);
        
        ylabel('Adjusted Network Strength Score', 'FontWeight', 'bold', 'FontSize', 14);
        
        % Main title + Subtitle with FDR stats
        title({clean_title, sprintf('(p_{FDR} = %.3f)', p_fdr(i))}, 'FontSize', 14, 'FontWeight', 'bold');
        
        % --- E. SAVE ---
        file_safe_name = strrep(clean_title, ' ', '_');
        file_safe_name = strrep(file_safe_name, '(', '');
        file_safe_name = strrep(file_safe_name, ')', '');
        
        out_name = fullfile(plot_dir, sprintf('%s_Violin_Plot.png', file_safe_name));
        exportgraphics(fig, out_name, 'Resolution', 300);
        close(fig); 
        
        fprintf('Saved significant plot: %s\n', out_name);
    end
    
    fprintf('\nAll significant plots successfully saved to the %s folder!\n', plot_dir);
end