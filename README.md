# EEG Signal Processing CHW2

## Overview

This project involves analyzing EEG signals with various processing techniques using MATLAB. The dataset includes multiple EEG files, and the tasks focus on signal visualization, dimensionality reduction, noise reduction, and artifact removal.

## Tasks

### 1. Three-Channel EEG Signal

- **Data**: `mat1.Ex` (EEG signal with 3 channels, sampling frequency 200 Hz)
  - **a)** Plot the three-channel EEG signal.
  - **b)** Use the `scatter3` function to visualize the 3D data and examine the data spread in different directions.
  - **c)** Apply PCA to the observations, compute the covariance matrix, and extract principal components and variances. Plot the principal directions and normalized data.
  - **d)** Repeat the PCA analysis using MATLAB’s PCA function.
  - **e)** Perform SVD on the observations and analyze the PCA-related concepts with the left and right singular matrices and singular values.

### 2. Interictal EEG Signal with Noise

- **Data**: `mat2.Ex` (32-channel interictal EEG signal, `org_X`), noisy signals (`1_noise_X` to `5_noise_X`)
  - **a)** Add noise with different SNRs (-10 dB and -20 dB) to the clean signal. Plot and compare the clean and noisy signals for each SNR.
  - **b)** Use PCA and ICA to extract sources. Compare results using methods such as `2Com`.
  - **c)** Identify and retain spike sources, discard others.
  - **d)** Reconstruct the denoised signal in the sensor domain and plot it.
  - **e)** Plot the denoised signals for channels 13 and 24, alongside the original and noisy signals.
  - **f)** Compute and visualize the relative RMSE (RRMSE) for each method and SNR.

### 3. EEG Signals with Artifacts

- **Data**: `1NewData` to `4NewData` (EEG signals from patients, sampling frequency 250 Hz)
  - **a)** Plot the signals in the time domain with channel labels using `disp_eeg.m` and `plotEEG.m`.
  - **b)** Examine noise and artifacts. Determine if ICA can effectively remove these artifacts.
  - **c)** Apply ICA to the signals and extract independent components and mixing matrix. Use `m.R2COM` for ICA.
  - **d)** Plot temporal, frequency, and spatial characteristics of each component. Use `m.pwelch` and `m.plottopomap` for frequency and spatial plots.
  - **e)** Store the indices of relevant sources in `SelSources` and reconstruct the denoised signal. Plot and compare it with the original signal. Assess if the source selection and denoising were effective.

## Requirements

- MATLAB
