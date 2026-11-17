%% Plot NNMF Connectivity Components in BrainNetViewer
close all; clear all; clc;

% Add paths to toolboxes
addpath('/path/to/toolboxes/brainnetviewer');
addpath('/path/to/toolboxes/brain-connectivity-toolbox/');

% Define input files for node and edge file creation
template_dir = '/path/to/template/files';
node_coords = load(fullfile(template_dir, 'atlas_coordinates.txt'));
node_labels_tbl = readtable(fullfile(template_dir, 'atlas_labels.txt'), 'Delimiter', '\t', 'ReadVariableNames', false);
node_labels = node_labels_tbl.Var1;
node_cat = load('/path/to/atlas_region_categories.txt');
options = '/path/to/brainnetoptions.mat';

% Base directory and updated path to your saved NNMF connectivity results
base_dir = '/path/to/project/directory/'; 
analysis_dir = fullfile(base_dir, 'Results_1', 'NNMF_Connectivity'); 
plot_out_dir = fullfile(base_dir, 'BrainNet_Plots');

if ~exist(plot_out_dir, 'dir')
    mkdir(plot_out_dir); 
end

num_rois = 38;

% Thresholds to render
thresholds = [0.05, 0.15];

% Define specific networks of interest per frequency band
target_networks = containers.Map();
target_networks('7-13Hz')  = [1];       % Alpha
target_networks('13-30Hz') = [1, 4];    % Beta
target_networks('30-55Hz') = [4, 6];    % Low Gamma
target_networks('65-80Hz') = [1];       % High Gamma

freq_keys = target_networks.keys;

for f = 1:length(freq_keys)
    freq_str = freq_keys{f};
    nets_to_plot = target_networks(freq_str);
    
    nnmf_file = fullfile(analysis_dir, sprintf('NNMF_Connectivity_Results_%s.mat', freq_str));
    
    if ~exist(nnmf_file, 'file')
        fprintf('File not found, skipping: %s\n', nnmf_file);
        continue; 
    end
    
    fprintf('\n==================================================\n');
    fprintf('Processing %s (Networks: %s)...\n', freq_str, mat2str(nets_to_plot));
    fprintf('==================================================\n');
    
    load(nnmf_file, 'W_final', 'optimal_k', 'upper_tri_idx');
    
    for c = nets_to_plot
        if c > optimal_k
            warning('Network %d requested for %s, but optimal_k is %d. Skipping.', c, freq_str, optimal_k);
            continue;
        end
        
        %% 1. RECONSTRUCT THE 38x38 ADJACENCY MATRIX
        conn_mat = zeros(num_rois, num_rois);
        conn_mat(upper_tri_idx) = W_final(:, c);
        conn_mat = conn_mat + conn_mat'; % Symmetrize
        
        %% 2. LOOP OVER SPECIFIED THRESHOLDS (5% and 15%)
        for th_idx = 1:length(thresholds)
            edge_threshold = thresholds(th_idx);
            th_tag = sprintf('thresh%02d', round(edge_threshold * 100)); % 'thresh05' or 'thresh15'
            
            % Proportional thresholding via BCT
            conn_mat_thresh = threshold_proportional(conn_mat, edge_threshold);
            
            %% 3. SAVE THE .EDGE FILE
            edge_file = fullfile(plot_out_dir, sprintf('Network_%d_%s_%s.edge', c, freq_str, th_tag));
            save(edge_file, '-ascii', 'conn_mat_thresh');
            
            %% 4. CREATE AND SAVE THE .NODE FILE
            bin_mat = double(conn_mat_thresh > 0);
            node_size = degrees_und(bin_mat)';
            node_size(node_size == 0) = 0.1; % Base size for disconnected nodes
            
            t = table();
            t.Var1 = node_coords(:, 1);
            t.Var2 = node_coords(:, 2);
            t.Var3 = node_coords(:, 3);
            t.Var4 = reshape(node_cat, [], 1);
            t.Var5 = node_size;
            t.Var6 = reshape(node_labels, [], 1);
            
            node_file = fullfile(plot_out_dir, sprintf('Network_%d_%s_%s.node', c, freq_str, th_tag));
            writetable(t, node_file, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', false);
            
            %% 5. PLOT AND EXPORT HIGH-RES IMAGE
            out_img = fullfile(plot_out_dir, sprintf('Connectivity_Network_%d_%s_%s.png', c, freq_str, th_tag));
            fprintf('Rendering Network %d (%s) at %d%% threshold -> %s\n', ...
                c, freq_str, round(edge_threshold * 100), out_img);
            
            % Generate BrainNet figure
            plotInBrainNet(edge_file, node_file, 'options', options);
            
            % Allow OpenGL pipeline to complete rendering
            drawnow; 
            pause(2);
            
            % Capture and write 300 DPI figure
            fig = gcf;
            exportgraphics(fig, out_img, 'Resolution', 300);
            close(fig);
        end
    end
end

fprintf('\nSuccess! Target connectivity networks rendered at 5%% and 15%% to: %s\n', plot_out_dir);