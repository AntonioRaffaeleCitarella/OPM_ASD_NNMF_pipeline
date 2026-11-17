# OPM_ASD_NNMF_pipeline

This repository contains the analytical pipeline used to process and analyze optically pumped magnetometer-based magnetoencephalography (OPM-MEG) data in young children (2–5 years old) with Autism Spectrum Disorder (ASD) and Typically Developing (TD) peers.
The pipeline implements a data-driven methodology to estimate regional cortical activity and functional connectivity, applies dimensionality reduction to extract shared network motifs, and runs robust statistical models to link these neural signatures to diagnostic and behavioral profiles.
  
Pipeline Overview.\
The codebase is organized into three main processing steps:
  1. Activity & Connectivity Estimation\
    1a. Activity: Computes the analytic signal via Hilbert transform and extracts the instantaneous amplitude envelope to estimate voxel-wise frequency-band activity (normalized via Coefficient of Variation) -> s7_opm_activity_downsampled-m\
    1b. Connectivity (AEC): Applies symmetric orthogonalization to regional time series to correct for spatial leakage, followed by Amplitude Envelope Correlation to estimate functional connectivity between regions\
  2. Dimensionality Reduction (PCA & NNMF) -> s9_pca_nnmf_activity.m and s_9_pca_nnmf_connectivity.m\
    2a. Data Preparation: Flattens activity maps and connectivity matrices, applying non-negativity thresholds.\
    2b. PCA / Truncated SVD: Reduces the feature space to retain 95% of the variance, preventing overfitting in subsequent steps.\
    2c. NNMF (Non-Negative Matrix Factorization): Extracts group-level spatial network topographies (Feature matrix, W) and participant-specific network strengths (Coefficient matrix, H) across different frequency bands.\
  3. Statistical Analysis -> General Linear Models to evaluate the relationship between NNMF component strengths (the H matrix) and:\
    3a. Diagnosis -> GLM_groupdifferences.m\
    3b. Non Verbal Developmental Quotient -> GLM_nvdq.m\
    3c. CBCL DSM-5 derived measures (subgroup ASD analysis) -> GLM_subgroup\

Requirements\
MATLAB
