function s7_opm_activity(p, sub_num, clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SCRIPT TO MULTIPLY WEIGHTS X TIMESERIES AND CALC ACTIVITY FOR NNMF
% % Run the study configuration
% p = opm_study_config_babyeinsteins_NNMF();
% Config based on the work and preprocesing pipeline by: Dr. Marlee M Vandewouw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the trial configuration json
if strcmp(p.epoching.type, 'json')
    fid = fopen(p.epoching.json);
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    trial_info = jsondecode(str);
end
%%
    % Only run if doesn't exist or if clobbering
    if ~exist([p.directories.sourcemodel_pre(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) '_activity_map_' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_MNI.nii.gz'], 'file') || clobber == 1
    
    % Now that all the channels have been identified, load in the bad
    % channels and perform MFC - do this for both the data and the noise
    fid = fopen(p.badchannels.bad_channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'r');
    raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    lines = strtrim(raw{1});
    lines = lines(~cellfun(@isempty, lines));
    
    % Split each line by comma, then flatten
    if isempty(lines)
        bad_chs = {};
    else
        split_lines = cellfun(@(x) strsplit(x, ','), lines, 'UniformOutput', false);
        bad_chs = strtrim(cat(1, split_lines{:}));
        bad_chs = bad_chs(~cellfun(@isempty, bad_chs));
    end

    % Read in the header
    hdr = ft_read_header(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'headerformat', 'opm_fil_mjt');
    
    % Load in the MEG data, CONTINUOUS, and demean
    cfg_loadds = [];
    cfg_loadds.dataset = p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num));
    cfg_loadds.headerformat = 'opm_fil_mjt';
    cfg_loadds.dataformat = 'opm_fil_mjt';
    cfg_loadds.eventformat = 'opm_fil_mjt';
    cfg_loadds.channel = hdr.grad.label;
    cfg_loadds.demean = 'yes';
    cfg_loadds.continuous = 'yes';
    
    % Load in the raw data
    data_continuous = ft_preprocessing(cfg_loadds);
    
    % BAD CHANNELS
    
    % Eliminate the bad channels
    cfg_ch = [];
    cfg_ch.channel = setdiff(data_continuous.label, bad_chs);
    data_continuous = ft_selectdata(cfg_ch, data_continuous);
    
    % Also apply to grad
    [~, ind] = intersect(data_continuous.grad.label, data_continuous.label, 'stable');
    grad_orig = data_continuous.grad.label;
    f = fieldnames(data_continuous.grad);
    for ii = 1:length(f)
        if length(data_continuous.grad.(sprintf('%s', f{ii}))) == length(grad_orig)
            data_continuous.grad.(sprintf('%s', f{ii})) = data_continuous.grad.(sprintf('%s', f{ii}))(ind, :);
        end
    end

    % INTERFERENCE SUPRESSION

    % Mean field correct
    N = data_continuous.grad.chanori;
    M = eye(length(N)) - N*pinv(N);
    data_continuous.trial{1} = M*data_continuous.trial{1};
        

    % FILTERING
    
    % Prepare preprocessing for bandpass filtering
    f = fieldnames(p.preproc);
    f = f(startsWith(f, 'bp'));
    cfg_preproc = [];
    for kk = 1:length(f)
        cfg_preproc.(sprintf('%s', f{kk})) = p.preproc.(sprintf('%s', f{kk}));
    end
    
    % Check if bandstop, and determine accordingly
    if isfield(p.preproc, 'bsfilter')
        if strcmp(p.preproc.bsfilter, 'yes')
            
            % Add to filter
            cfg_preproc.bsfilter = 'yes';
            
            if isnumeric(p.preproc.bsfreq)
                
                % Add the bandstop
                cfg_preproc.bsfreq = [p.preproc.bsfreq'-p.preproc.bswidth, p.preproc.bsfreq'+p.preproc.bswidth];
                
            elseif strcmp(p.preproc.bsfreq, 'noise_peaks')
                
                % Load in the peak widths
                cfg_preproc.bsfreq = readmatrix(p.filter.bandstop_widths(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
                
            end
        end
    end
    
    % Now filer
    try
        data_continuous_filt = ft_preprocessing(cfg_preproc, data_continuous);
    catch
        cfg_preproc.bsfiltord = 3;
        data_continuous_filt = ft_preprocessing(cfg_preproc, data_continuous);
    end
    
    % ICA
    
    % Reject bad ICA components
    cfg = [];
    if exist(p.ica.bad_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')
        % Load bad components
        cfg.component = load(p.ica.bad_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        % Load ICA components
        components = load(p.ica.ica_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        % Reject
        data_continuous_filt = ft_rejectcomponent(cfg, components, data_continuous_filt);
    end
    
    % Store
    data_continuous_filt_orig = data_continuous_filt;
    data_continuous_filt.orig = data_continuous_filt_orig;
    
    % RE-FILTER
    
    % Re-filter: grid
    if contains(p.sourcemodel.method, 'grid')
        f = fieldnames(p.preproc);
        f = f(startsWith(f, 'bp'));
        cfg_preproc2 = [];
        for kk = 1:length(f)
            cfg_preproc2.(sprintf('%s', f{kk})) = p.preproc.(sprintf('%s', f{kk}));
        end
        cfg_preproc2.bpfreq = p.sourcemodel.grid.foi;
        data_continuous_filt.grid = ft_preprocessing(cfg_preproc2, data_continuous_filt.orig);
    end
    
    % EPOCH
    
    % Load in the trial configuration
    cfg_trl = load(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
    
    % Epoch
    cfg_epoch = [];
    cfg_epoch.trl = cfg_trl.trl;
    data_filt.orig = ft_redefinetrial(cfg_epoch, data_continuous_filt.orig);

    data_filt.grid = ft_redefinetrial(cfg_epoch, data_continuous_filt.grid);
    
    % PSEUDOCONTINUOUS

    sensor_data = cat(2,data_filt.grid.trial{:});


    % LOAD IN WEIGHTS

    load(p.sourcemodel.grid.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num),p.task(sub_num)),'W_grid');

    % MULTIPLY WEIGHTS x DATA

    grid_ve = W_grid' * sensor_data;

    % COMPUTE HILBERT ENVELOPE
    
    % Note: MATLAB's hilbert() function operates on columns by default. 
    % Because grid_ve is [voxels x time], we must transpose it ('), apply hilbert, 
    % and transpose it back so it doesn't mix up voxels with time points.
    analytic_sig = hilbert(grid_ve')';
    
    % Extract the magnitude to get the amplitude envelope
    grid_env = abs(analytic_sig);

    % GET ACTIVITY MEASURE (STD/MEAN)
    
    mu = mean(grid_env, 2);
    sigma = std(grid_env, 0, 2);
    activity_cv = sigma ./ mu;

    % CONVERT INTO NIFTI

    % Load
    meshes = load([p.headmodel.template_dir(sub_num) '/meshes.mat']);
    
    % Extract
    mri_segmented = meshes.segmentedmri;
    
    % Create the headmodel
    cfg_hm = [];
    cfg_hm.method = p.headmodel.method;
    vol = ft_prepare_headmodel(cfg_hm, mri_segmented);

    % Convert to the source position units
    vol = ft_convert_units(vol, p.sourcemodel.units);

    % Prepare the source model - in template space
    cfg_sm = [];
    cfg_sm.grid.resolution = p.sourcemodel.grid.res;
    cfg_sm.grid.unit = p.sourcemodel.units;
    cfg_sm.grid.tight = 'yes';
    cfg_sm.headmodel = vol;
    sourcemodel_grid = ft_prepare_sourcemodel(cfg_sm);
 
    % Convert to the same units
    sourcemodel_grid = ft_convert_units(sourcemodel_grid, vol.unit);

    % Put in FieldTrip structure in order to write nifti
    pow = [];
    pow.Q = activity_cv;
    pow.pos = sourcemodel_grid.pos;
    pow.unit = sourcemodel_grid.unit;
    pow.dim = sourcemodel_grid.dim;
 
    % Write (Fixed double semicolon and added -f to gzip)
    cfg_wn = [];
    cfg_wn.filename =  [p.directories.sourcemodel_pre(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) '_activity_map_' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz.nii'];
    cfg_wn.parameter = 'Q';
    cfg_wn.precision = 'double';
    cfg_wn.filetype = 'nifti';
    ft_sourcewrite(cfg_wn, pow);
    cmd = ['gzip -f ' p.directories.sourcemodel_pre(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) '_activity_map_' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz.nii'];
    system(cmd);

    % TRANSFORM INTO MNI SPACE
    setenv('LD_LIBRARY_PATH', strrep(getenv('LD_LIBRARY_PATH'), '/usr/local/MATLAB/R2024a/sys/os/glnxa64:', ''));

    % Register to MNI space
    cmd = ['antsApplyTransforms -d 3 -i ' p.directories.sourcemodel_pre(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) '_activity_map_' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz.nii.gz -r ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_Warped.nii.gz -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_1Warp.nii.gz -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_0GenericAffine.mat -o ' p.directories.sourcemodel_pre(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) '_activity_map_' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_MNI.nii.gz'];
    system(cmd);
    end

    % =========================================================
    % NORMALIZATION & DOWNSAMPLING BLOCK
    % =========================================================
    out_basename = [p.directories.sourcemodel_pre(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) '_activity_map_' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_MNI_norm'];
    
    if ~exist([out_basename '.nii.gz'], 'file') || clobber == 1
        
        % 1. Load the high-res 1mm MNI map generated by ANTs
        mni_activity_map = ft_read_mri([p.directories.sourcemodel_pre(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) '_activity_map_' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_MNI.nii.gz']);

        % 2. DOWNSAMPLE TO 6MM (Matches your true MEG grid resolution)
        fprintf('Downsampling MNI map to 6mm...\n');
        cfg = [];
        cfg.resolution = 6; 
        mni_activity_map = ft_volumereslice(cfg, mni_activity_map);

        % 3. Normalize (Doherty method)
        activity_map = mni_activity_map.anatomy;
        activity_map_norm = correct_global_effects(activity_map);
        mni_activity_map.anatomy = activity_map_norm;

        % 4. Write as .nii first, THEN gzip it (prevents FieldTrip crashes)
        ft_write_mri([out_basename '.nii'], mni_activity_map, 'dataformat', 'nifti');
        system(['gzip -f ' out_basename '.nii']);
        
    end

end

function [adjusted_data] = correct_global_effects(data)
    % Extracts all valid, non-zero values (ignores ANTs tiny interpolation noise)
    vals = data(~isnan(data) & abs(data) > 1e-6 & ~isinf(data));
    
    % 1. Find the "noise peak" (the mode of the distribution)
    [f, xi] = ksdensity(vals);
    [~, max_idx] = max(f);
    peak_mu = xi(max_idx);
    
    % 2. Estimate initial variance robustly (Median Absolute Deviation)
    robust_std = 1.4826 * median(abs(vals - peak_mu));
    
    % 3. Isolate the data within +/- 1 SD of the peak
    core_noise_data = vals(vals >= (peak_mu - robust_std) & vals <= (peak_mu + robust_std));
    
    % 4. Fit a Gaussian to this specific noise peak to get the final null parameters
    pd = fitdist(core_noise_data, 'Normal');
    null_mu = pd.mu;
    null_sigma = pd.sigma;
    
    % 5. Z-score the original data using the null noise distribution
    adjusted_data = (data - null_mu) ./ null_sigma;
    
    % 6. MAKE STRICTLY NON-NEGATIVE (Zero-clipping for NNMF)
    adjusted_data(adjusted_data < 0) = 0;
    
    % Restore the exact zeros (e.g., empty space outside the brain)
    adjusted_data(data == 0) = 0;
    adjusted_data(isnan(data)) = NaN;

end