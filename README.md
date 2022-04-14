# iRD-paper

### MATLAB Code for the iRD Paper, submitted to NeuroImage 2022

This code belongs to the preprint "Indirect decoding of behavior from brain signals via reconstruction of experimental conditions" by Soch & Haynes (2022), publicly available from *bioRxiv* and submitted to *NeuroImage*. It consists of MATLAB scripts processing fMRI data, as available from OpenNeuro, into results and figures, as available from the paper.

- Preprint: https://www.biorxiv.org/content/10.1101/2022.03.24.485588
- Data: https://openneuro.org/datasets/ds001734 (Version 1.0.4)
- Code: https://github.com/JoramSoch/iRD-paper


## Requirements

This code was developed and run using the following software:
- Windows 10
- [MATLAB R2020a](https://de.mathworks.com/help/matlab/release-notes.html)
- [MATLAB Statistics Toolbox](https://de.mathworks.com/help/matlab/release-notes.html)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (Revision 7771 as of 14/01/2020)
- [ITEM Toolbox](https://github.com/JoramSoch/ITEM) (as on GitHub)
- [The Decoding Toolbox](https://sites.google.com/site/tdtdecodingtoolbox/) (Version 3.999)
- [Bayesian Prevalence](https://github.com/robince/bayesian-prevalence) package (MATLAB)


## Instructions

For re-running analyses of the empirical data, you need to perform the following steps:
1. Create a folder on your computer that hereafter is referred to as the "study directory".
2. Download the empirical data set from https://openneuro.org/datasets/ds001734 and place it into a sub-folder of the study directory called "data".
3. Download the analysis scripts from https://github.com/JoramSoch/iRD-paper and place them into a sub-folder of the study directory called "tools".
4. Open MATLAB, set your current directory to this sub-folder "tools", edit the study directory in the `project_directories.m` and run this script.
5. Run the script `analysis_MGT_iRD.m` located in the same folder. Ideally, run the code step by step to ensure that each step of the analysis succeeds.
6. When the entire analysis of the empirical data set has been performed, the script `Figure_4_5_6.m` from the tools directory should generate Figures 4/5/6, as they appear in the paper.

The sub-folder "NBD_tools/" contains functions that are used by `analysis_MGT_iRD.m`, `MGT_decode_XZ_wb_TDT.m` and `MGT_show_results_wb_TDT.m` to perform the indirect response decoding (iRD).


## Bonus: Graphical Abstract

<img src="https://github.com/JoramSoch/iRD-paper/raw/master/Figure_GA.png" alt="Graphical Abstract" width=1000>

This graphical abstract illustrates the core idea of the paper: Behavioral responses can either be predicted from the experimental design (stimulus-based response decoding, sbRD), directly decoded from the stimulus-related fMRI signals (direct response decoding, dRD) or indirectly decoded by first reconstructing the experimental design from the stimulus-related fMRI signals (indirect response decoding, iRD). Click <a href="https://github.com/JoramSoch/iRD-paper/blob/master/Figure_GA.pdf">here</a> for a PDF version.