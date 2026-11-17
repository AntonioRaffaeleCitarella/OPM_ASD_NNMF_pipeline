close all; clear; clc;

% 1. Add BrainNet Viewer to path
bnvRoot = '/path/to/toolboxes/brainnetviewer';
addpath(genpath(bnvRoot));

% 2. Define the Brain Surface Template
surf = fullfile(bnvRoot, 'Data', 'SurfTemplate', 'BrainMesh_ICBM152_smoothed.nv');

% 3. Define the NIFTI Activity File (Corrected Path)
base_dir = '/path/to/project/directory/';
nifti_gz = fullfile(base_dir, 'activity_nifti_file.nii.gz');
nifti_unzipped = fullfile(base_dir, 'pca_nifti_file.nii');

% Check if the file is zipped, and safely unzip it for BrainNet Viewer
if exist(nifti_gz, 'file') && ~exist(nifti_unzipped, 'file')
    fprintf('Unzipping NIfTI file for BrainNet Viewer...\n');
    gunzip(nifti_gz);
end

% 4. Create an Output Directory for Activity Plots (Optional, currently commented out)
plot_out_dir = '/path/to/output/directory';
if ~exist(plot_out_dir, 'dir')
    mkdir(plot_out_dir); 
end
out_img = fullfile(plot_out_dir, 'Output_Activity_Network_1.png');

% 5. Execute BrainNet Viewer
fprintf('Loading Network 2 into BrainNet Viewer...\n');

% Pass the UNZIPPED .nii file to the viewer
BrainNet_MapCfg(surf, nifti_unzipped);