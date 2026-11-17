function s9_pca_nnmf_activity(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s9_pca_nnmf_activity(p)
%
% 1) Aggregates 3D NIfTI activity maps across all participants, 
% 2) Applies a group brain mask, 
% 3) Runs Truncated-SVD Denoising (95% variance threshold)
% 4) Runs NNMF to identify group-level spatial networks 
% 5) Reconstructs the networks back into 3D NIfTI volumes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_subs = size(p.subject_data, 1);

% Define the path to your master spreadsheet
excel_path = '/path/to/nnmf_results.xlsx';

% Loop through every frequency band defined in the config file
for ff_idx = 1:size(p.sourcemodel.grid.foi, 1)
    
    % Get the current frequency limits for naming
    freq_band = p.sourcemodel.grid.foi(ff_idx, :);
    freq_str = sprintf('%d_%dHz', freq_band(1), freq_band(2));

    freq_prefix = sprintf('x%d_%d', freq_band(1), freq_band(2));
    
    fprintf('\n=========================================================\n');
    fprintf('Running Activity NNMF for Frequency Band: %s\n', freq_str);
    fprintf('=========================================================\n');

    %% 1. AGGREGATE NIFTI ACTIVITY MAPS AND APPLY BRAIN MASK
    fprintf('Loading participant Activity NIfTIs...\n');

    % Find the first valid subject file to extract the 3D dimensions & header
    first_sub_file = '';
    for ss = 1:num_subs
        temp_file = [p.directories.sourcemodel_pre(p.subject(ss), p.session(ss), p.run(ss), p.task(ss)) '_activity_map_' freq_str '_MNI_norm.nii.gz'];
        if exist(temp_file, 'file')
            first_sub_file = temp_file;
            break;
        end
    end
    
    if isempty(first_sub_file)
        warning('No NIfTI files found for %s! Skipping to next frequency.', freq_str);
        continue;
    end
    
    % Read the header info
    nii_info = niftiinfo(first_sub_file);
    vol_dim = nii_info.ImageSize; 

    % Preallocate a temporary 4D matrix to hold everyone's 3D brain
    % Dimensions: [X, Y, Z, Participants]
    all_vols = zeros([vol_dim, num_subs]);

    valid_subs = 0;
    for ss = 1:num_subs
        % Define file path
        nii_file = [p.directories.sourcemodel_pre(p.subject(ss), p.session(ss), p.run(ss), p.task(ss)) '_activity_map_' freq_str '_MNI_norm.nii.gz'];
        
        if exist(nii_file, 'file')
            % Load the 3D volume
            vol = niftiread(nii_file);
            
            % Ensure no negatives (just in case) and replace NaNs with 0
            vol(vol < 0) = 0; 
            vol(isnan(vol)) = 0;
            
            all_vols(:, :, :, ss) = vol;
            valid_subs = valid_subs + 1;
        else
            warning('Activity NIfTI missing for Subject %d.', ss);
        end
    end

    % CREATE THE BRAIN MASK
    % Find all voxels that have a value > 0 in at least one participant
    mean_vol = mean(all_vols, 4);
    brain_mask = mean_vol > 0;
    num_voxels = sum(brain_mask(:)); 

    fprintf('Masked out empty space. Found %d valid brain voxels.\n', num_voxels);

    % Flatten the 3D brains into a 2D matrix for NNMF [Voxels x Participants]
    V = zeros(num_voxels, num_subs);
    for ss = 1:num_subs
        current_vol = all_vols(:, :, :, ss);
        V(:, ss) = current_vol(brain_mask); % Extracts only the valid voxels into a column
    end

    clear all_vols % Free up computer memory!
    
    %% 1.5 TRUNCATED SVD DENOISING 
    fprintf('\nStarting Truncated-SVD Denoising for 3D Activity Maps...\n');
    
    % Run SVD (Economy size to save memory)
    [U_svd, S_svd, V_svd] = svd(V, 'econ');
    
    % Calculate variance explained by each component
    eigenvalues = diag(S_svd).^2;
    variance_explained = eigenvalues / sum(eigenvalues);
    cumulative_variance = cumsum(variance_explained);
    
    % Find the cutoff for 95% variance
    variance_threshold = 0.95; 
    num_pca_comps = find(cumulative_variance >= variance_threshold, 1);
    
    fprintf('Keeping top %d spatial components (explains %.1f%% of variance).\n', ...
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

    % Define a tiny threshold to determine if a component is "represented" 
    rep_threshold = 1e-4; 

    fprintf('\nStarting iterative NNMF to find optimal k on denoised data...\n');

    for k = 2:max_k

        % Run NNMF
        [W, H] = nnmf(V, k, 'Replicates', 10, 'Options', statset('MaxIter', 1000));


        % Create a binary matrix: 1 if represented, 0 if not
        is_represented = H > rep_threshold;

        % Calculate the percentage of participants representing EACH component
        prop_per_comp = sum(is_represented, 2) / num_subs;

        % RULE 1: Each component must be in at least 50% of participants
        rule1_passed = all(prop_per_comp >= 0.50);

        % RULE 2: Mean representation across all components must be >= 70%
        rule2_passed = mean(prop_per_comp) >= 0.70;

        % Print progress
        fprintf('Testing k=%d | Min Rep: %.1f%% | Mean Rep: %.1f%% ... ', ...
            k, min(prop_per_comp)*100, mean(prop_per_comp)*100);

        % Evaluate rules
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


    %% 4. RECONSTRUCT AND SAVE THE NETWORKS AS 3D NIFTIS

    output_dir = fullfile(p.directories.base_dir, 'analysis', 'NNMF_Activity');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    fprintf('\nSaving Group Network NIfTIs...\n');

    for c = 1:optimal_k

        % 1. Create a blank 3D volume matching header datatype
        network_3D = zeros(vol_dim, nii_info.Datatype);

        % 2. Push the 1D network voxels exactly back into the brain mask space
        network_3D(brain_mask) = cast(W_final(:, c), nii_info.Datatype);

        % 3. Save as NIfTI using the header info from the first subject
        out_filename = fullfile(output_dir, sprintf('Group_Activity_Network_70_ica_pca_%d_%s.nii', c, freq_str));
        niftiwrite(network_3D, out_filename, nii_info);

        % Compress it
        system(['gzip -f ' out_filename]);

        fprintf('Saved Network %d to %s\n', c, output_dir);

    end

    % Save the H_matrix, W_matrix, optimal k, and the rebuilding index
    save(fullfile(output_dir, sprintf('NNMF_Activity_Results_70_ica_pca_%s.mat', freq_str)), 'H_final', 'W_final', 'optimal_k', 'brain_mask');
    fprintf('\nDone! Frequency %s pipeline complete.\n', freq_str);

    %% 5. EXPORT RESULTS TO SPREADSHEET
    if exist(excel_path, 'file')
        fprintf('Appending results to %s...\n', excel_path);

        % Load existing spreadsheet
        results_table = readtable(excel_path);

        % Append columns for each component found
        for c = 1:optimal_k
            % Create column names matching your previous format
            col_name_raw = sprintf('%s_Activity_%d_Raw_70_ica_pca', freq_prefix, c);

            % Insert data into the table
            results_table.(col_name_raw) = H_final(c, :)';
        end

        % Save updated spreadsheet
        writetable(results_table, excel_path);
        fprintf('Successfully added %d new columns to spreadsheet.\n', optimal_k);
    else
        warning('Could not find %s! Make sure the Excel file exists.', excel_path);
    end

end % End of frequency loop

fprintf('\nSuccess! All Activity NNMF processing is finished.\n');

end