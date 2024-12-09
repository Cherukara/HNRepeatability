#!/bin/bash

# This script will loop through the HN data and create transformation matrices
# that align each scan with the first session (based on the magnitude data), 
# and then apply this to the QSMs
#
# Matthew Cherukara, October 2024

# Directory
datadir=/media/cherukara/DATA/HN_Repeatability_BIDS

# Method name
method=VSHARP

# Loop through subjects
for s in {1..10}
do
    # Define subject number
    if [ $s -lt 10 ]
    then
        ss=0${s}
    else
        ss=${s}
    fi

    # Define more specific directories
    anatdir=${datadir}/rawdata/sub-${ss}/ses-01/anat
    qsmdir=${datadir}/derivatives/qsm/sub-${ss}/ses-01/qsm

    # Define subject name for session 1
    oname="sub-${ss}_ses-01"

    # # Create a brain-only QSM for session 1
    #fslmaths ${qsmdir}/${oname}_unwrapped-SEGUE_bfr-PDF_susc-autoNDI_Chimap \
    #    -mul ${anatdir}/${oname}_desc-brain_mask \
    #    ${qsmdir}/${oname}_mask-brain_susc-${method}_Chimap

    # Loop over sessions
    for i in {2..6}
    do
        # Define the name of the scan
        sname="sub-${ss}_ses-0${i}"

        # Print the stage we are at
        echo "Processing $sname"

        # Define the directories
        sanatdir=${datadir}/rawdata/sub-${ss}/ses-0${i}/anat
        sqsmdir=${datadir}/derivatives/qsm/sub-${ss}/ses-0${i}/qsm

        # Create a brain-only QSM for this session
        fslmaths ${sqsmdir}/${sname}_unwrapped-SEGUE_mask-nfe_bfr-VSHARP_susc-autoNDI_Chimap \
            -mul ${sanatdir}/${sname}_desc-brain_mask \
            ${sqsmdir}/${sname}_mask-brain_method-${method}_Chimap

        # # Register this brain to session 1 (rigid-body only)
        #flirt -in ${sqsmdir}/${sname}_mask-brain_susc-${method}_Chimap \
        #    -ref ${qsmdir}/${oname}_mask-brain_susc-autoNDI_Chimap \
        #    -out ${sqsmdir}/${sname}_mask-brain_reg-ses01_susc-${method}_Chimap \
        #    -omat ${sqsmdir}/regbrain_${i}to1.mat \
        #    -dof 6

        # Apply already-existing transformation matrix
        flirt -in ${sqsmdir}/${sname}_mask-brain_method-${method}_Chimap \
            -ref ${qsmdir}/${oname}_mask-brain_susc-autoNDI_Chimap \
            -out ${sqsmdir}/${sname}_mask-brain_reg-ses01_method-${method}_Chimap \
            -init ${sqsmdir}/regbrain_${i}to1.mat \
            -applyxfm

        # Re-apply the original mask
        fslmaths ${sqsmdir}/${sname}_mask-brain_reg-ses01_method-${method}_Chimap \
            -mul ${anatdir}/${oname}_desc-brain_mask \
            ${sqsmdir}/${sname}_mask-brain_reg-ses01_method-${method}_Chimap

    done

done
