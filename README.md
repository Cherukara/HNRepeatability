# HNRepeatability

Copyright Matthew Cherukara, 9 December 2024. See LICENSE.

MATLAB code accompanying the paper
- Cherukara MT, Shmueli K. Comparing repeatability metrics for quantitative susceptibility mapping in the head and neck. *Magn. Reson. Mater. Phy.* 2025 [doi.org/10.1007/s10334-025-01229-3](https://doi.org/10.1007/s10334-025-01229-3)


## Dependencies

- STI-Suite 3.0 [Berkeley](https://people.eecs.berkeley.edu/~chunlei.liu/software.html)  
- MEDI-Toolbox [Cornell](http://pre.weill.cornell.edu/mri/pages/qsm.html)  
- FANSI-Toolbox [CMilovic](https://gitlab.com/cmilovic/FANSI-toolbox)  
- chi-separation toolbox [SNU-LIST](https://github.com/SNU-LIST/chi-separation)

You will need to edit **qsm_pipeline.m** with paths to your local installation of STI-Suite and MEDI-Toolbox. FANSI and chi-separation are linked as submodules to this one. The **utils/** folder contains other necessary functions. Please properly credit the authors of any toolboxes used.

## How To Use

Requires multi-echo complex data, stored as NIFTIs, in BIDS format, with separate files for each echo and for the magntiude and phase data. Such as:
```
dataset
└── rawdata 
    ├── sub-01
        │   ├── ses-01
        │   │   └── anat
        │   │       ├── sub-01_ses-01_part-mag_echo-1_GRE.nii
        │   │       ├── sub-01_ses-01_part-mag_echo-2_GRE.nii
        │   │       ├── sub-01_ses-01_part-phase_echo-1_GRE.nii
        │   │       ├── sub-01_ses-01_part-phase_echo-2_GRE.nii
        │   │       └── etc.
        │   └── ses-02
        │       └── anat
        │           └── etc.
        └── sub-02
            └── etc.
```

First, use **qsm_pipeline.m** to reconstruct QSMs using the methods specified within that script. The outputs will be saveed in a BIDS-formatted **derivatives/** folder within your BIDS directory.

For repeatability analysis, you will also need a brain mask for each subject, and a ROI segmentation map. 

For voxel-wise repeatability metrics (NRMSE and XSIM), you will need to register the QSMs of subsequent repetitions to the first session. This is done using the bash script **align_hn.sh**.

Once the QSMs are aligned, run **voxel_repeatability.m** to calculate the repeatability metrics in each ROI for each subject. The results will be stored as a **.mat** file. You can then use **analysis_voxelwise.m** to display the results (reproducing Figure 3 from the paper).

For ROI-based repeatability metrics (standard deviations, ICC, etc.), run **roi_averages.m**, which will load your reconstructed QSMs and save out the average susceptibility values for each ROI into a **.mat** file.

Use **analysis_repeatability.m** and **analysis_anova.m** to calculate repeatability metrics on the resulting data. These scripts can reproduce Figures 4-8 in the paper.

## Head and Neck QSM Repeatability Data

This code was designed to be used with the Head and Neck QSM Repeatability Data set, which is available at
- Cherukara MT, Shmueli K. Head and Neck QSM Repeatability Data. *University College London*: Dataset. 2024. [doi.org/10.5522/04/27993215](https://doi.org/10.5522/04/27993215)