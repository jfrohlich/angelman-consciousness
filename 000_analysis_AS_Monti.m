% 000_analysis
% Joel Frohlich
% Last updated: 17 November, 2022
 
000_analysis_AS_Monti.m % This file describes each script and the analysis pipeline

%% General overview of work flow

# For the new paper (machine learning), data are first imported using 
# ana_preproc_as_*_butterHP_firLP.m, where "*" designates the specific data
# being imported (e.g., "Dup15q_wsleep"). Imported data are manually 
# inspected and ICA is performed using AS_data_scroll.m, and data are then
# postprocessed using ana_AS_postproc_butterhp_firlp.m,
# ana_AS_postproc_Dup15q.m, or ana_AS_postproc_TD.m. After running the
# appropriate analysis scripts to compute EEG features, one should run 
# ana_AS_2021_ML_and_STATS to train on AS data and validate on healthy and
# dup15q data ... this script will then do stats and figures. Finally, for
# the NT replication (train on healthy data), one should run 
# ana_AS_2021_ML_Train_On_TD.m

## directories

./2021_analysis            # Output for the new paper (machine learning)
./2021_analysis/TrainOnTD  # Output specific to the neurotypical replication (train on healthy data) in the new paper
./scripts                  # These are the scripts and functions described in this file
./scripts/archived         # Old scripts, consider deleting 

## wrappers and helpers 
AS_data_scroll.m               # calls scroll_data.m, contains code for AS specific options, such as labels for delta patterns 
plot_topo.m                    # wrapper function for fieldtrip plot_topo
plot_topo_AS.m                 # AS specific wrapper for topo plotting 
plot_topo_AS_classic.m         #  This wrapper function calls the "classic" (i.e., original) FT topoplot function with thin countour lines 
effect_size.m                  # helper function: cohen's d 
makefigpretty.m                # Renders figure in format suitable for journal publication (This has been moved to Universal directory)
ana_AS_find_sections_PHZ.m     # Finds appropraite sleep/wake sections to match using the findchangepts() function and taking into account posterior hot zone (PHZ) in sleep 
ana_AS_sleep_PHZ               # This is the function that tests our method of using PHZ to classify sleep into dreaming and no expereince. Use this to generate figures for reviewers.
ana_AS_max_pow_wake.m          # Identifies sections of awake EEG to match with sleep on power (by finding the awake EEG with highest delta power)
myresample.m                   # Fixes bug in the native Matlab function
myround.m                      # Used for finding smoothing windows with gMLZ; this is a wrapper version of the native MATLAB round function that allows user to specify whether we should round towards nearest ODD integer. 
gpu2cpu.m                      # gathers data from GPU arrays back into regular arrays so they can be read on computer without usable GPU
powernoise.m                   # Generates power-law noise (not used by any analysis scripts but useful for playing with new functions) 
Double_check_annotaions.m      # Double check that we don't have sleep marked incorrectly (e.g., that there are no awake tags in the sections marked as sleep)
freq_wavelet_bias_test.m       # Test to see if allow_nan parameter biases the frequency transform 
sleep_field_repair.m           # Honestly don't remember exactly ... something to do with check the sleep fields
MSE_debugging                  # Test effects of NaNs on the ro_mse tool
LZC_debugging                  # Test effects of NaNs on the ro_LZC tool 
ana_AS_update_fields_butterHPfirLP_2021.m # This take the already filtered data (FIR lowpass, IIR (butterworth) highpass) and populate the datastructure fields with info from the FIR bandpass filtered data, including ICA weights.
CalcPermEn.m                   # Helper function to compute the permutation entropy, calls code by Jaco Sitt
cohen_d.m                      # Computes effect sizes as Cohen's d
compare_annotations.m          # compares our sleep annotations with those done by a trained neurologist through Ben Philpot's lab
Compare_filter_responses.m     # Compare the performance of FIR and Butterworth filters, used to choose a filter
Compare_filters.m              # Slightly more streamlined than the above script, this one loads data that have already been filtered and plots the powerspectra
FDRCorrectAndCombineTables.m   # For the Tables in the main manuscript and supplement reporting p-values for AUCs (either comparisons of entropy vs spectral or better than chance performance), this will do the FDR correction on all p-values and then generate the tables. 
FDRCorrectAndCombineTabesTrainOnTD.m # Same as above, but for the replication on healthy data
Find_drowsy_tags.m             # After loading a dataset, this will look for drowsy tags and related annotations (e.g., yawning)
get_annotations.m              # Computes the degree of overlap in sleep annotations (our annotations vs Philpot annotations)
Read_Nespeca_Annotations.m     # Reads the annotations by Mark Nespeca and writes them to a CSV file
ReportMeds.m                   # This is the script that generated the supplemental table of medications for the first paper in Neurosci of Consci
Transfer_sleep_labels          # Move the sleep labels from old ICA files to the files with new preprocessing/filtering
Update_Nespeca_annotations.m   # This adds sleep blocks from a table to the existing data field:  data.cfg.dattype.AWAKE_NESPECA = [data.cfg.dattype.AWAKE_NESPECA; T.awakeStart(irow)*data.fsample T.awakeStop(irow)*data.fsample];
ro_dwPLI.m                     # Test the debiased weighted phase lag index code on some toy data


## Frequency transforms
ro_freq_wavelet_JF.m       # This is the main frequency transform. The "_JF" is added at the end so that the original file created by Joerg wouldn't accidentally be overwritten
ro_freq_wavelet_TFT.m      # This is the freq transforms that does time-frequency analysis (no averaging across time). Used for finding sections of high wakeful power
freq_wavelet_bias_test.m   # Quantify the bias induced by allowing NaNs in the frequency transform window 

## Tools

ro_LZC.m                   # Computes the gMLZ Lempel-Ziv from Yeh 2018, as well as the "vanilla" Lempel-Ziv
ro_LZCv.m                  # Only compute the vanilla Lempel-Ziv, all other fields will be left empty. This is the one used for the new paper (machine learning)
ro_mse.m                   # Computes MSE with parameters to choose from Xie or Costa, can set dyanamic tolerance
ro_PermEn.m                # Computes the permutation entropy by calling code by Jaco Sitt's group
ro_wPLI.m                  # Compute the debiased weighted phase lag index
ro_wSMI.m                  # Computes the weighted symbolic information by calling code by Sitt's group
surrogate_data.m           # Crates surrogate data with same amplitude distribution using FFT phase randomization (calls code that D Toker gave me)

## Preprocessing

ana_preproc_as_boston_butterHP_firLP.m # Applies 5th order Butterworth IIR filter as highpass and FIR filter a lowpass to boston files for AS (discontinuous EDFs, excludes EDF C)
ana_preproc_as_boston__EDF_C_butterHP_firLP.m # Applies 5th order Butterworth IIR filter as highpass and FIR filter a lowpass to boston AS files (EDF C)
ana_preproc_as_Dup15q_wsleep_butterHP_firLP.m # Applies 5th order Butterworth IIR filter as highpass and FIR filter a lowpass to Dup15q files
ana_preproc_AS_SD_2019_butterHP_firLP.m # Applies 5th order Butterworth IIR filter as highpass and FIR filter a lowpass to San Diego AS files
ana_preproc_AS_TD_wsleep_butterHP_firLP.m # Applies 5th order Butterworth IIR filter as highpass and FIR filter a lowpass to healthy TD controls


## postprocessing
ana_AS_postproc_butterhp_firlp.m   # This is the current postproc script that removes ICA artifacts, interpolates channels, excludes drowsy sections. It follows a preproc script that applies an FIR lowpass filter and an IIR (Butterorth) highpass filter that minimally attenuates the 0.5 - 1.0 Hz slow band
ana_AS_posproc_Dup15q.m            # The current postproc script for Dup15q (remove ICA artifacts, interpolate bad channels, etc.)
ana_AS_postproc_TD.m               # The current postproc script for healthy controlles (remove ICA artifacts, interpolate bad channels, etc.)

## analysis 

ana_AS_2021_ML_and_STATS.m         # This is the main script that runs the stats and generates figures for the newest paper (Machine learning)
ana_AS_2021_ML_and_STATS_Jeff.m    # This is a slight modification of the main analysis script that implements a suggestion by Jeff for the stats (currently unused, move to deletion candidates?)
ana_AS_2021_ML_Train_On_TD.m       # This is just like ana_AS_2021_ML_and_STATS, EXCEPT that it implements the neurotypical replication analysis (i.e., the training is done on healthy children rather than Angelman syndrome)
ana_AS_criticality.m               # This is the script to compute D. Toker's criticality measure on AS data. Unused in our paper because criticality is not an entropy measure--our analysis focused on spectral vs entropy
ana_AS_criticality_Dup15q.m        # Computes Toker's criticality measure on Dup15q data--not used for reasons mentioned for ana_AS_criticality.m
ana_AS_criticality_TD.m            # Computes Toker's criticality measure on healthy data--not used for reason mentioned above.
ana_AS_FT_allow_nan_0dot2.m        # Computes frequency transform for AS with allow_nan = 0.2 THIS IS THE VERSION WE USED FOR BOTH PAPERS
ana_AS_FT_Dup15q_allow_nan_0dot2.m # Same as above for Dup15q
ana_AS_FT_TD_allow_nan_0dot2.m     # Same as above for healthy controls
ana_AS_LZC.m                       # Computes Lempel-Ziv for AS
ana_AS_LZC_Dup15q.m                # Same as above for Dup15q
ana_AS_LZC_TD_controls.m           # Same as above for healthy controls
ana_AS_MSE_Xie_dynr.m              # Computes multiscale entropy for AS with dynamic tolerance (we used this version)
ana_AS_MSE_Xie_dynr_Dup15q_controls.m # Same as above for Dup15q
ana_AS_MSE_Xie_dynr_TD_controls.m  # Same as above for healthy controls
ana_AS_PermEn_Decomp.m             # This calls the permutation entropy decomposition on the AS data and does the stats, figures, etc.
ana_AS_stats_cluster.m             # Performed stats (including cluster permutation statistics) and generated figures for the first paper (Neurosci of Consci)
ana_AS_stats_cluster_short.m       # Abbreviated version of the above script. For all stats/figures, run the full script above.
ana_AS_TFT.m                       # This performs a time-frequency transform on the data, e.g., for finding sections of high/low delta power in the first paper (Neurosci of Consci)
ana_AS_wPLI_fieldtrip.m            # Compute the dibiased weighted phase locking index on AS data using Fieldtrip
ana_AS_wPLI_TD_fieldtrip.m         # Same as above, but for healthy controls
ana_AS_wPLI_Dup15q.m               # Same as above, but for Dup15q
ana_AS_wSMI.m                      # Routine for the weighted symbolic mutual information on AS data by calling code by Jaco Sitt
ana_AS_wSMI_TD.m                   # Same as above, but for healthy controls
ana_AS_wSMI_Dup15q.m               # Same as above, but for Dup15q 
CountPeaksFooof.m                  # Count spectral peaks (channel-averaged) for each participant
MSE_simulated_data.m               # Does MSE on toy data--useful for understanding MSE on different data types
PermEnDecomposition.m              # Modified code from Pedro, this is what actually run the Permutation Entropy decomposition 

## Stats and plotting

AS_big_picture_figure.m                # Plots barplots of ROC AUCs for each feature type and group
connectivity_distance_plotting.m       # Generates sup. figure 7 for the new figure, showing histogra of channel distances and matrix of channel-pair groupings
Dup15q_example_EEG_plotting.m          # Plots Dup15q EEG traces for Fig. 1 in the new paper
Figure_1_New.m                         # Generates subpanels for Fig. 1 in the first paper (Neurosci of Consci)
LZc_illustration.m                     # Generates schematic of LZC for Fig. 1 in the new paper (machine learning)
SimulatedEEG.m                         # Simulated EEG signals for the conceptual figure in the discussion section of AS paper (Neurosci of Consci 2020)
Supplemental_EEG_Figure                # Figure requested by reviewer showing longer sections of EEG and power spectra for sleep wake in 27-month-old girl with AS
permclustertest.m                      # Does permutation cluster statistics on UNPAIRED data
RMpermclusttest.m                      # Does permutation cluster statistics on paired or repeated measures data
ana_AS_stats_cluster.m                 # Does stats and figures using permutation cluster statistics
ana_AS_stats_cluster_short.m           # Same as ana_AS_stats_cluster.m, but with old/unneeded sections omitted (just the figures that should go into paper)
LZC_power_intrasubject_correlation.m   # Does the correlations within subject (30 s windows) between LZC and delta power
MSE_power_intrasubject_correlation.m   # Does the correlations within subject (30 s windows) between MSE and delta power


## Deletion candidates (these are now archived or deleted)

ana_AS_ANOVA.R              # Unused R script that looks at each EEG measure with ANOVA (IVs: sleep, group, age)
ana_AS_CD.m                 # Compute the causal density (Anil Seth)
ana_AS_CD_mo.m              # Find the proper model order for causal density
ana_AS_CD_mo_awake.m        # Find the proper model order (awake EEG only)
ana_AS_coh.m                # Unused script, precessor to the dwPLI script ... began writing this one to look at coherence between channels
ana_AS_DE.m                 # Differential entropy for AS, this measure didn't work well on our data
ana_AS_DE_TD.m              # Differential entropy for healthy controls, this measure didn't work well on our data
ana_AS_DFA.m                # Compute DFA (unused)
ana_AS_FT.m                 # Does the frequency transform calling ro_freq_wavelet_JF, allow_nan = 0.5
ana_AS_LZC_decomp.m         # Abandoned attempt at Lempel-Ziv decomposition. The newest paper instead does a similar decomp on permutation entropy 
ana_AS_MSE_Xie.m            # Computes the multiscale entropy for AS -- this version is archive because it does NOT dynamicaly adjust the tolerance r for each timescale
ana_AS_PCA.m                # This was for a PCA of the AS EEG --Does PCA on the awake and asleep data concatenated together
ana_AS_phi.m                # Compute phi from EEG
ana_AS_phi_sleep.m          # Compute phi from sleep EEG (already deleted?)
ana_AS_phi_awake.m          # Compute phi from awake EEG (already deleted?)
ana_AS_postproc.m           # Data postprocessing: remove bad ICs, interpolate bad channels, reref to average
ana_AS_stats_controls.m            # An early version of the stats script for the new paper--this was replaced by ana_AS_2021_ML_and_STATS.m 
ana_AS_stats_controls_2021.m       # Another early version of the stats script for the new paper -- replaced by ana_AS_2021_ML_and_STATS.m 
ana_AS_update_fields_butter_2021.m # Populated the data sturcture fields with info from the old, FIR filtered data.
eegfordipfit.m              # reformat the data so that we can do dipfit (dipole source localization). This was never used.

ana_AS_DFA_stats.m        # stats for DFA -- hasn't be used
ana_AS_LZC_stats.m        # was the stats for lempel ziv, not used anymore
ana_AS_match_power.m      # old method of matching wake/sleep power 
ana_AS_MSE.m              # Does MSE on all AS data by calling ro_mse.m--does not specify dynamic tolerance
correlation_phenotype_unused.m         # Shows correlations between gMLZ and all covaries (unused in manuscript)
ro_freq_wavelet_dfa.m     # Unused -- this was written for the DFA analysis 
ro_DE.m                   # Implements the differential entropy, never used
progressbar.m             # Downloaded from fileexchange, never used 
test_regression_models.m  # Not sure if I remember what this one was for? 





