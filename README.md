
# 2024 Dynamic Mode Brain fMRI

This repository accompanies the manuscript titled "[Large-scale signal propagation modes in the human brain](https://www.sciencedirect.com/science/article/pii/S105381192500360X)."

## Dynamic Modes (DMs)
<p align="center">
  <strong>Principal DM (DM1)</strong> (unimodal–transmodal axis) &nbsp;&nbsp; | &nbsp;&nbsp; 
  <strong>Bi-asym DM (DM5)</strong> (inter-hemispheric propagation)
</p>

<div align="center">
  <img src="./visualize/Supplementary_Video_1.gif" 
       alt="Principal DM - DM1"
       title="Principal DM (DM1): propagation along unimodal–transmodal axis"
       width="350" />
  <img src="./visualize/Supplementary_Video_5.gif" 
       alt="Bi-asym DM - DM5"
       title="Bi-asym DM (DM5): inter-hemispheric propagation"
       width="350" />
</div>

## Prerequisites

### MATLAB Toolboxes

-   **SPM12** – [website](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
    
-   **conn** – [website](https://web.conn-toolbox.org/)
    
-   **cifti-matlab** – [GitHub Repository](https://github.com/Washington-University/cifti-matlab.git)
    
-   **lag-code** – [GitHub Repository](https://github.com/ryraut/lag-code)
    
-   **Accelerated Permutation Inference for the ACE Model** – [GitHub Repository](https://github.com/NISOx-BDI/APACE)
    
-   **Analysis code for the paper "Exploring the Latent Structure of Behavior Using the Human Connectome Project’s Data"** – [GitHub Repository](https://github.com/connectomicslab/hcp-behavioral-domains)  
    _Note: To properly run this code, the numbers in the file names within the processing scripts should be removed._
    
-   **DataViz** – For violin plots. [GitHub Repository](https://github.com/povilaskarvelis/DataViz)
    

### Atlases

-   **Glasser 2016** – [Download Link](https://balsa.wustl.edu/study/RVVG)
    
-   **The Cole-Anticevic Brain-wide Network Partition** – [GitHub Repository](https://github.com/ColeLab/ColeAnticevicNetPartition)


### Human Connectome Project (HCP) Data

The following HCP S1200 data are required for running the analyses in this repository:

- **Resting-state fMRI**
  - ICA-FIX denoised resting-state fMRI data (REST1 & REST2)
  - Minimally preprocessed resting-state fMRI data (REST1 & REST2)

> [!NOTE]
> For reproducibility, it may be helpful to know exactly which datasets were used in our analyses:
>   - Resting-state fMRI 1 FIX-Denoised (Extended)  
>   - Resting-state fMRI 2 FIX-Denoised (Extended)  
>   - Resting-state fMRI 1 Minimally Preprocessed  
>   - Resting-state fMRI 2 Minimally Preprocessed  
>   *The Compact version of the FIX-Denoised dataset should also work, since only files included in the Compact release were used. However, the folder structure may differ, so minor modifications to the scripts might be required.*  

- **Task-state fMRI**
  - Minimally preprocessed task fMRI data (7 tasks: WM, EMOTION, MOTOR, LANGUAGE, GAMBLING, SOCIAL, RELATIONAL)  

- **Demographic and additional data**
  - Behavioral data (`.csv` file, unrestricted dataset)  
  - FreeSurfer data (`.csv` file, unrestricted dataset)  
  - Restricted data (`.csv` file, includes genetic information)  

Access to restricted and genetic data requires approval from the HCP database.
    

## File Descriptions

-   **Code00_Setup.m** : (**Requires modification!**) Setup script for saving dataset paths that must be downloaded from HCP. *To properly run the other scripts, you first need to execute this setup script after modifying the dataset path variables according to your computing environment.*

-   **Code01_REST_*** : Scripts for analyzing resting-state fMRI data
    
    -   **Code01_REST_01_Extracting_ROI_timeseries.m** – Extract ROI time series from fMRI data in CIFTI format
        
    -   **Code01_REST_02_Filtering.m** – Bandpass filtering of the extracted ROI time series
        
    -   **Code01_REST_03_DM_analysis_ext_fbDMD_CV.m** – Cross-validation for predicting future BOLD signals using varying numbers of DMs, previous ML benchmarks (linear model & manifold-based model), and a null model
        
    -   **Code01_REST_04_DM_test_retest_reliability.m** – Test-retest reliability assessment using HCP REST1 and REST2 data (results not included in the paper)
        
    -   **Code01_REST_05_DM_analysis_ext_fbDMD.m** – Group-level DMD using the entire dataset and estimation of DM metrics for all subjects. The algorithm is identical to *Code01_REST_03*
        
    -   **Code01_REST_06_rsFMRI_properties.m** – Derivation of canonical resting-state features in fMRI data
        
    -   **Code01_REST_07_dFC_analysis.m** – Time-varying FC analysis
        
    -   **Code01_REST_08_DM_specification.m** – Evaluation of each DM by calculating spatial correlations with canonical features and testing reconstruction ability
        
    -   **Code01_REST_09_Heritability_analysis.m** – Heritability assessment using the ACE model
        
    -   **Code01_REST_10_Behavior_analysis_all.m** – Correlation analysis between DM metrics and latent behavioral factors estimated by [Schöttner et al. (2023)](https://www.nature.com/articles/s41598-022-27101-1)
        
    -   **Code01_REST_11_DM_display_all.m** – Generate animated videos of DM evolution
        
    -   **Code01_REST_12_Null_Test_Presence_DM.m** – Test for the presence of DMs in given fMRI data
        
-   **Code02_Task_*** : Scripts for analyzing task-state fMRI data
    
    -   **Code02_Task_01_Extracting_ROI_timeseries_tfMRI.m** – Extract ROI time series from fMRI data in CIFTI format
        
    -   **Code02_Task_02_Extracting_Task_Info.m** – Extract onset and duration of task execution
        
    -   **Code02_Task_03_Filtering_CompCor_tfMRI.m** – Denoising of task-state data with bandpass filtering, regression of movement regressors, detrending, and event-related activation
        
    -   **Code02_Task_04_DM_analysis_indiv_fitting.m** – Estimation of DM metrics from task-state data using group-level DM obtained from resting-state dataset
        
    -   **Code02_Task_05_DM_similarity.m** – Assessment of DM metric similarity across tasks or subjects
        
    -   **Code02_Task_06_DM_difference.m** – Assessment of DM metric differences between tasks
        
    -   **Code02_Task_07_DM_significance.m** – Test for the presence of DMs in given fMRI data
        
    -   **Code02_Task_08_Summarizing_Subject_Info.m** – Script for organizing demographic information used in this study
        
