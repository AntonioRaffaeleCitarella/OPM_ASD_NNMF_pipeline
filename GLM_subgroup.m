% 1. LOAD CONFIGURATION AND SPREADSHEET
p = opm_study_config_babyeinsteins_NNMF();

% Pointing to the new Results_1 folder where the updated spreadsheet lives
table_file = '/path/to/nnmf_results.xlsx'; 

% Load master table
data_master = readtable(table_file);

% Format diagnosis and SUBSET TO ASD ONLY
data_master.dx = categorical(strtrim(data_master.dx));
data_asd = data_master(data_master.dx == 'ASD', :);

fprintf('Filtered data to ASD subgroup only (n = %d total ASD participants).\n\n', height(data_asd));

% Define the target CBCL traits to loop through
target_traits = {'ADOS'}; 
%target_traits = {'anxp', 'adhdp', 'autsp', 'oddp', 'dp'}; 

% Loop through each trait
for t = 1:length(target_traits)
    current_trait = target_traits{t};
    fprintf('========================================================================================\n');
    fprintf('STARTING ASD-ONLY ANALYSIS FOR: %s\n', upper(current_trait));
    fprintf('========================================================================================\n');
    
    % Start fresh with the ASD subset for this trait
    data = data_asd;
    
    % Ensure trait is numeric and DROP MISSING DATA for this specific score
    if ~isnumeric(data.(current_trait))
        data.(current_trait) = str2double(data.(current_trait));
    end
    valid_idx = ~isnan(data.(current_trait));
    data = data(valid_idx, :);
    
    % PREPARE COVARIATES
    data.sex = categorical(strtrim(data.sex)); 
    
    % Ensure age is numeric and Z-SCORE IT (based on this specific valid subset)
    if ~isnumeric(data.age)
        data.age = str2double(data.age);
    end
    data.age = (data.age - mean(data.age, 'omitnan')) / std(data.age, 'omitnan');
    
    fprintf('Running Doherty Robust GLM on %d ASD participants with valid %s scores...\n', height(data), upper(current_trait));
    
    % IDENTIFY ALL RAW ICA COLUMNS 
    all_vars = data.Properties.VariableNames;
    raw_cols = all_vars(contains(all_vars, '_Raw_70_ica_pca')); 
    
    num_components = length(raw_cols);
    if num_components == 0
        error('No raw components found in the table. Check column names!');
    end
    
    % INITIALIZE STORAGE ARRAYS
    std_beta_trait = zeros(num_components, 1);
    ci_lower       = zeros(num_components, 1);
    ci_upper       = zeros(num_components, 1);
    p_vals         = zeros(num_components, 1);
    outliers_drop  = zeros(num_components, 1);
    p_age          = zeros(num_components, 1);
    p_sex          = zeros(num_components, 1);
    models         = cell(num_components, 1);
    
    % RUN ROBUST LINEAR MODELS
    for i = 1:num_components
        col_name = raw_cols{i};
        
        % Z-SCORE THE NETWORK STRENGTH
        y_raw = data.(col_name);
        data.y_std = (y_raw - mean(y_raw, 'omitnan')) / std(y_raw, 'omitnan');
        
        % Define the formula dynamically (dx is removed since everyone is ASD)
        formula = sprintf('y_std ~ %s + sex + age', current_trait);
        
        % INITIAL FIT
        mdl_initial = fitlm(data, formula);
        
        % COOK's DISTANCE OUTLIER REJECTION
        cooksD = mdl_initial.Diagnostics.CooksDistance;
        mean_cooks = mean(cooksD, 'omitnan');

        outlier_idx = (cooksD > (3 * mean_cooks)) | (abs(cooksD) > 0.5);
        outliers_drop(i) = sum(outlier_idx);

        data_clean = data;
        data_clean(outlier_idx, :) = [];
        
        % FINAL ROBUST FIT
        mdl_final = fitlm(data_clean, formula);
        models{i} = mdl_final;
        
        % Extract Beta and p-value for the current CBCL trait
        std_beta_trait(i) = mdl_final.Coefficients{current_trait, 'Estimate'};
        p_vals(i)         = mdl_final.Coefficients{current_trait, 'pValue'};
        
        % Extract Confidence Intervals
        ci = coefCI(mdl_final);
        trait_idx = strcmp(mdl_final.CoefficientNames, current_trait);
        ci_lower(i) = ci(trait_idx, 1);
        ci_upper(i) = ci(trait_idx, 2);
        
        % Extract Covariate p-values
        p_age(i) = mdl_final.Coefficients{'age', 'pValue'};
        
        sex_coef_name = mdl_final.CoefficientNames(contains(mdl_final.CoefficientNames, 'sex_'));
        if ~isempty(sex_coef_name)
            p_sex(i) = mdl_final.Coefficients{sex_coef_name{1}, 'pValue'};
        else
            p_sex(i) = NaN; 
        end
    end
    
    % SEPARATE FDR (BENJAMINI-HOCHBERG) CORRECTION 
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
    
    % DISPLAY FINAL SUMMARY TABLE
    col_names = {'Network', sprintf('Std_Beta_%s', upper(current_trait)), 'CI_95_Low', 'CI_95_High', ...
                 sprintf('p_%s_Raw', upper(current_trait)), 'p_FDR', 'FDR_Sig', ...
                 'Tests_Corrected_For', 'Outliers_Removed', 'p_Age', 'p_Sex'};
             
    summary_table = table(raw_cols', std_beta_trait, ci_lower, ci_upper, p_vals, p_fdr, fdr_sig, tests_applied, outliers_drop, p_age, p_sex, ...
        'VariableNames', col_names);
    
    disp(summary_table);
    
    % Save results
    report_path = fullfile('/path/to/results/output/', sprintf('Stats_Robust_GLM_ASD_Only_%s_Raw_FDR.csv', upper(current_trait)));
    writetable(summary_table, report_path);
    fprintf('\nSaved %s report to: %s\n', upper(current_trait), report_path);
    
    % 8. GENERATE PUBLICATION-READY SCATTER PLOTS (SIGNIFICANT ONLY)
    fprintf('\nGenerating Publication-Ready Scatter Plots for SIGNIFICANT results...\n');
    plot_dir = fullfile('/path/to/results/output/', sprintf('GLM_Publication_Plots_ASD_Only_%s_Raw_FDR', upper(current_trait)));
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
            
            if contains(raw_name, 'Connectivity')
                mod_label = 'Connectivity Network';
                tokens = regexp(raw_name, 'Connectivity_(\d+)', 'tokens');
            else
                mod_label = 'Activity Component';
                tokens = regexp(raw_name, 'Activity_(\d+)', 'tokens');
            end
            
            if ~isempty(tokens)
                net_num = tokens{1}{1};
            else
                net_num = 'X';
            end
            
            clean_title = sprintf('%s %s %s', freq_label, mod_label, net_num);
            
            % --- B. EXTRACT PARTIAL RESIDUALS (Adjusted Scores) ---
            b0 = mdl.Coefficients{'(Intercept)', 'Estimate'};
            b_trait = mdl.Coefficients{current_trait, 'Estimate'};
            x_val = mdl.Variables.(current_trait);
            
            % The adjusted Y mathematically strips away age and sex
            y_adj = mdl.Residuals.Raw + b0 + (b_trait .* x_val);
            
            % --- C. DRAW THE CUSTOM SCATTER PLOT ---
            fig = figure('Name', clean_title, 'Color', 'w', 'Position', [200, 200, 600, 500], 'Visible', 'off');
            hold on;
            
            % Orange for ASD
            color_ASD = [0.8500 0.3250 0.0980]; 
            
            % Plot the Scatter Dots
            scatter(x_val, y_adj, 70, color_ASD, 'filled', 'MarkerFaceAlpha', 0.8, 'HandleVisibility', 'off');
            
            % Plot the Main Effect Regression Line
            x_limits = [min(x_val), max(x_val)];
            y_line = b0 + b_trait .* x_limits;
            plot(x_limits, y_line, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            
            % --- D. FORMATTING & AESTHETICS ---
            set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off');
            
            ylabel('Adjusted Network Strength Score', 'FontWeight', 'bold', 'FontSize', 13);
            xlabel(sprintf('%s Score', upper(current_trait)), 'FontWeight', 'bold', 'FontSize', 13);
            
            % Main title + Subtitle with FDR p-value
            title({clean_title, sprintf('(p_{FDR} = %.3f)', p_fdr(i))}, 'FontSize', 14, 'FontWeight', 'bold');
            
            % Clean up the legend
            legend('off');
            
            % --- E. SAVE ---
            file_safe_name = strrep(clean_title, ' ', '_');
            file_safe_name = strrep(file_safe_name, '(', '');
            file_safe_name = strrep(file_safe_name, ')', '');
            
            out_name = fullfile(plot_dir, sprintf('%s_%s_Scatter.png', file_safe_name, upper(current_trait)));
            exportgraphics(fig, out_name, 'Resolution', 300);
            close(fig); 
            
            fprintf('Saved significant plot: %s\n', out_name);
        end
    end
    fprintf('Completed all plots for %s.\n\n', upper(current_trait));
end

fprintf('SUCCESS! All 5 ASD-Only CBCL trait analyses are fully complete.\n');