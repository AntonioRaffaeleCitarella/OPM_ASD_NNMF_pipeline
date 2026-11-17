function s8_pca_nnmf_connectivity(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s8_pca_nnmf_connectivity(p)
% 
% Aggregates 38x38 connectivity matrices across all participants, applies 
% SVD denoising, and runs Non-Negative Matrix Factorization (NNMF) to 
% identify group-level networks.
% 
% Updates: 
% - Loads _sprem (median spike removal) connectivity matrices.
% - INTEGRATED TRUNCATED SVD DENOISING (90% Variance Threshold)
% - Appends H_final results automatically to NNMF_results.xlsx

% % Run the study configuration
% p = opm_study_config_babyeinsteins_NNMF();
% Config based on the work and preprocesing pipeline by: Dr. Marlee M Vandewouw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_subs = size(p.subject_data, 1);
num_rois = 38; % Based on  atlas

% Define the path to your master spreadsheet
excel_path = '/path/to/master/spreadsheet.xlsx';

% Loop through every frequency band you defined in the config file
for ff_idx = 1:size(p.connectivity.fois, 1)
    
    % Get the current frequency limits for naming
    freq_band = p.connectivity.fois(ff_idx, :);
    freq_str = sprintf('%d-%dHz', freq_band(1), freq_band(2));
    
    % Create a clean prefix for MATLAB table columns (e.g., 'x4_7')
    freq_prefix = sprintf('x%d_%d', freq_band(1), freq_band(2));
    
    fprintf('\n=========================================================\n');
    fprintf('Running NNMF for Frequency Band: %s\n', freq_str);
    fprintf('=========================================================\n');

    %% 0. LOAD THE DATA
    all_adjmats = zeros(num_rois, num_rois, num_subs);
    valid_subs = 0;
    
    for ss = 1:num_subs
        % Dynamically point to the right frequency index (ff_idx)
        orig_file = p.conenctivity.aec_mat_ve_coordinates(p.subject(ss), p.session(ss), p.run(ss), p.task(ss), 'all', sprintf('%d', ff_idx));
        
        % Alter the filename to load the _sprem version
        [fpath, fname, fext] = fileparts(orig_file);
        mat_file = fullfile(fpath, [fname '_sprem' fext]);
        
        if exist(mat_file, 'file')
            temp_load = load(mat_file, 'adjmat');
            all_adjmats(:, :, ss) = temp_load.adjmat;
            valid_subs = valid_subs + 1;
        else
            warning('Connectivity matrix missing for Subject %d. Check s7b!', ss);
        end
    end

    if valid_subs == 0
        warning('No data found for %s. Skipping to next frequency.', freq_str);
        continue;
    end
    
    fprintf('Successfully stacked %d participants.\n', valid_subs);

    %% 1. PREPARE THE DATA FOR NNMF

    % Find the indices for the upper triangle (703 unique edges)
    upper_tri_idx = find(triu(ones(num_rois, num_rois), 1));
    num_edges = length(upper_tri_idx);

    % Initialize the V matrix [703 edges x N participants]
    V = zeros(num_edges, num_subs);

    % Flatten each participant's matrix into a column of V
    for ss = 1:num_subs
        current_mat = all_adjmats(:, :, ss);
        V(:, ss) = current_mat(upper_tri_idx);
    end

    % BIOLOGICAL CLIP: NNMF strictly requires non-negative data.
    % This removes anti-correlations so SVD focuses only on positive networks.
    V(V < 0) = 0; 
    
    %% 1.5 TRUNCATED SVD DENOISING 
    fprintf('\nStarting Truncated-SVD Denoising...\n');
    
    % Run SVD (Economy size to save memory)
    [U_svd, S_svd, V_svd] = svd(V, 'econ');
    
    % Calculate variance explained by each component
    eigenvalues = diag(S_svd).^2;
    variance_explained = eigenvalues / sum(eigenvalues);
    cumulative_variance = cumsum(variance_explained);
    
    % Find the cutoff for 95% variance
    variance_threshold = 0.95; 
    num_pca_comps = find(cumulative_variance >= variance_threshold, 1);
    
    fprintf('Keeping top %d components (explains %.1f%% of variance).\n', ...
        num_pca_comps, cumulative_variance(num_pca_comps) * 100);
    
    % Reconstruct the clean matrix using ONLY the top components
    U_trunc = U_svd(:, 1:num_pca_comps);
    S_trunc = S_svd(1:num_pca_comps, 1:num_pca_comps);
    V_trunc = V_svd(:, 1:num_pca_comps);
    
    V_denoised = U_trunc * S_trunc * V_trunc';
    
    % MATHEMATICAL MICRO-CLIP: 
    % Fixes tiny negative artifacts (e.g., -0.00001) caused by orthogonal rotation.
    num_neg = sum(V_denoised(:) < 0);
    fprintf('Clipping %d microscopic negative reconstruction ripples to 0...\n', num_neg);
    V_denoised(V_denoised < 0) = 0;
    
    % Replace the raw V matrix with our clean, denoised V matrix
    V = V_denoised;

    %% 2. ITERATIVELY FIND OPTIMAL k COMPONENTS

    max_k = 20; % Set a reasonable upper limit to test
    optimal_k = 2; % Start at 2 components
    rep_threshold = 1e-4; 

    fprintf('\nStarting iterative NNMF to find optimal k on denoised data...\n');

    for k = 2:max_k
        % Run NNMF (using replicates to avoid local minima)
        [W, H] = nnmf(V, k, 'Replicates', 10, 'Options', statset('MaxIter', 1000));
        
        
        % Create a binary matrix: 1 if represented, 0 if not
        is_represented = H > rep_threshold;
        
        % Calculate the percentage of participants representing EACH component
        prop_per_comp = sum(is_represented, 2) / num_subs;
        
        % EVALUATE DOHERTY RULES
        rule1_passed = all(prop_per_comp >= 0.50);
        rule2_passed = mean(prop_per_comp) >= 0.70; 
        
        fprintf('Testing k=%d | Min Rep: %.1f%% | Mean Rep: %.1f%% ... ', ...
            k, min(prop_per_comp)*100, mean(prop_per_comp)*100);
        
        if rule1_passed && rule2_passed
            fprintf('PASSED\n');
            optimal_k = k; % Update our best k
        else
            fprintf('FAILED\n');
            fprintf('Stopping criteria met. Optimal components: k = %d\n', optimal_k);
            break; % Exit the loop
        end
    end

    %% 3. RUN FINAL NNMF WITH OPTIMAL k

    fprintf('\nRunning final NNMF reconstruction with k = %d...\n', optimal_k);
    [W_final, H_final] = nnmf(V, optimal_k, 'Replicates', 50, 'Options', statset('MaxIter', 2000));

    % Save the .mat results securely using the config paths
    out_folder = p.directories.base_dir; 
    final_filename = fullfile(out_folder, sprintf('NNMF_Connectivity_Results_%s_sprem_70.mat', freq_str));

    save(final_filename, 'H_final', 'W_final', 'optimal_k', 'upper_tri_idx');
    fprintf('Saved .mat results to: %s\n', final_filename);
    
    %% 4. EXPORT RESULTS TO SPREADSHEET
    if exist(excel_path, 'file')
        fprintf('Appending results to %s...\n', excel_path);
        
        % Load existing spreadsheet
        results_table = readtable(excel_path);
        
        
        % Append columns for each component found
        for c = 1:optimal_k
            % Create column names matching your previous format
            col_name_raw = sprintf('%s_Connectivity_%d_sprem_Raw_70_ica_pca', freq_prefix, c);
            
            % Insert data into the table
            results_table.(col_name_raw) = H_final(c, :)';
        end
        
        % Save updated spreadsheet
        writetable(results_table, excel_path);
        fprintf('Successfully added %d new columns to spreadsheet.\n', optimal_k * 2);
    else
        warning('Could not find %s! Make sure the Excel file exists.', excel_path);
    end
    
end % Ends the frequency loop

fprintf('\nDone! Group-level NNMF pipeline complete for all frequencies.\n');
end