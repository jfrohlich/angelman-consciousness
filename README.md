000_analysis
Joel Frohlich
Last updated: 17 November, 2022
 
### Preamble ###
This was a big project with lots of scripts. Confused? If you want to replicate 
the machine learning analysis, the script you're looking for is ana_AS_2021_ML_and_STATS.m, 
which runs the actual stats and machine learning on EEG features. If you are 
looking to replicate other parts of the project, or just wanting to understand
the dependencies of the main script, see below for a list of what each .m file
actually does. 

## General overview of work flow

* Data are imported using ana_preproc_as_* _butterHP_firLP.m, where "* " designates the specific data
being imported (e.g., "Dup15q_wsleep"). 
* Imported data are manually inspected and ICA is performed using AS_data_scroll.m, and data are then
postprocessed using ana_AS_postproc_butterhp_firlp.m, ana_AS_postproc_Dup15q.m, or ana_AS_postproc_TD.m. 
* After running the appropriate analysis scripts to compute EEG features, one should run ana_AS_2021_ML_and_STATS to train on AS data and validate on healthy and dup15q data ... this script will then do stats and figures. Finally, for the NT replication (train on healthy data), one should run 
ana_AS_2021_ML_Train_On_TD.m

## directories

* ./2021_analysis            # Output for the new paper (machine learning)
* ./2021_analysis/TrainOnTD  # Output specific to the neurotypical replication (train on healthy data) in the new paper
* ./scripts                  # These are the scripts and functions described in this file
* ./scripts/archived         # Old scripts, consider deleting 

## wrappers and helpers 
*  plot_topo.m                    # wrapper function for fieldtrip plot_topo
*  plot_topo_AS.m                 # AS specific wrapper for topo plotting 
*  plot_topo_AS_classic.m         #  This wrapper function calls the "classic" (i.e., original) FT topoplot function with thin countour lines 
*  makefigpretty.m                # Renders figure in format suitable for journal publication (This has been moved to Universal directory)
*  myresample.m                   # Fixes bug in the native Matlab function
*  myround.m                      # Used for finding smoothing windows with gMLZ; this is a wrapper version of the native MATLAB round function that allows user to specify whether we should round towards nearest ODD integer. 
*  gpu2cpu.m                      # gathers data from GPU arrays back into regular arrays so they can be read on computer without usable GPU
*  ana_AS_update_fields_butterHPfirLP_2021.m # This take the already filtered data (FIR lowpass, IIR (butterworth) highpass) and populate the datastructure fields with info from the FIR bandpass filtered data, including ICA weights.
*  CalcPermEn.m                   # Helper function to compute the permutation entropy, calls code by Jaco Sitt
*  cohen_d.m                      # Computes effect sizes as Cohen's d

## Tools

* ro_LZC.m                   # Computes the gMLZ Lempel-Ziv from Yeh 2018, as well as the "vanilla" Lempel-Ziv
* ro_LZCv.m                  # Only compute the vanilla Lempel-Ziv, all other fields will be left empty. This is the one used for the new paper (machine learning)
* ro_mse.m                   # Computes MSE with parameters to choose from Xie or Costa, can set dyanamic tolerance
* ro_PermEn.m                # Computes the permutation entropy by calling code by Jaco Sitt's group
* ro_wPLI.m                  # Compute the debiased weighted phase lag index
* ro_wSMI.m                  # Computes the weighted symbolic information by calling code by Sitt's group
* surrogate_data.m           # Crates surrogate data with same amplitude distribution using FFT phase randomization (calls code that D Toker gave me)
* ro_freq_wavelet_JF.m       # This is the main frequency transform. The "_JF" is added at the end so that the original file created by Joerg wouldn't accidentally be overwritten

## Preprocessing

* ana_preproc_as_boston_butterHP_firLP.m # Applies 5th order Butterworth IIR filter as highpass and FIR filter a lowpass to boston files for AS (discontinuous EDFs, excludes EDF C)
* ana_preproc_as_boston__EDF_C_butterHP_firLP.m # Applies 5th order Butterworth IIR filter as highpass and FIR filter a lowpass to boston AS files (EDF C)
* ana_preproc_as_Dup15q_wsleep_butterHP_firLP.m # Applies 5th order Butterworth IIR filter as highpass and FIR filter a lowpass to Dup15q files
* ana_preproc_AS_SD_2019_butterHP_firLP.m # Applies 5th order Butterworth IIR filter as highpass and FIR filter a lowpass to San Diego AS files
* ana_preproc_AS_TD_wsleep_butterHP_firLP.m # Applies 5th order Butterworth IIR filter as highpass and FIR filter a lowpass to healthy TD controls


## postprocessing
* ana_AS_postproc_butterhp_firlp.m   # This is the current postproc script that removes ICA artifacts, interpolates channels, excludes drowsy sections. It follows a preproc script that applies an FIR lowpass filter and an IIR (Butterorth) highpass filter that minimally attenuates the 0.5 - 1.0 Hz slow band
* ana_AS_posproc_Dup15q.m            # The current postproc script for Dup15q (remove ICA artifacts, interpolate bad channels, etc.)
* ana_AS_postproc_TD.m               # The current postproc script for healthy controlles (remove ICA artifacts, interpolate bad channels, etc.)

## analysis 

* ana_AS_2021_ML_and_STATS.m         # This is the main script that runs the stats and generates figures for the newest paper (Machine learning)
* ana_AS_2021_ML_Train_On_TD.m       # This is just like ana_AS_2021_ML_and_STATS, EXCEPT that it implements the neurotypical replication * analysis (i.e., the training is done on healthy children rather than Angelman syndrome)
* ana_AS_FT_allow_nan_0dot2.m        # Computes frequency transform for AS with allow_nan = 0.2 
* ana_AS_FT_Dup15q_allow_nan_0dot2.m # Same as above for Dup15q
* ana_AS_FT_TD_allow_nan_0dot2.m     # Same as above for healthy controls
* ana_AS_LZC.m                       # Computes Lempel-Ziv for AS
* ana_AS_LZC_Dup15q.m                # Same as above for Dup15q
* ana_AS_LZC_TD_controls.m           # Same as above for healthy controls
* ana_AS_MSE_Xie_dynr.m              # Computes multiscale entropy for AS with dynamic tolerance (we used this version)
* ana_AS_MSE_Xie_dynr_Dup15q_controls.m # Same as above for Dup15q
* ana_AS_MSE_Xie_dynr_TD_controls.m  # Same as above for healthy controls
* ana_AS_PermEn_Decomp.m             # This calls the permutation entropy decomposition on the AS data and does the stats, figures, etc.
* ana_AS_wPLI_fieldtrip.m            # Compute the dibiased weighted phase locking index on AS data using Fieldtrip
* ana_AS_wPLI_TD_fieldtrip.m         # Same as above, but for healthy controls
* ana_AS_wPLI_Dup15q.m               # Same as above, but for Dup15q
* ana_AS_wSMI.m                      # Routine for the weighted symbolic mutual information on AS data by calling code by Jaco Sitt
* ana_AS_wSMI_TD.m                   # Same as above, but for healthy controls
* ana_AS_wSMI_Dup15q.m               # Same as above, but for Dup15q 
* CountPeaksFooof.m                  # Count spectral peaks (channel-averaged) for each participant
* PermEnDecomposition.m              # Modified code from Pedro, this is what actually run the Permutation Entropy decomposition 
