% VOXEL_REPEATABILITY.m calculates voxel-wise repeatability metrics.
%   Created by MT Cherukara, October 2024
%
%
%       Copyright (C) University College London, 2024
%
%
% Takes reconstructed QSM data in a BIDS-compliant format, with specified ROIs,
% and calculates some voxel-wise measures of repeatability (NRMSE and XSIM).
% These results are then stored as a .mat file for later analysis. 
%
% QSM daa should have been reconstructed using QSM_PIPELINE.m, and then the
% results co-registered across repetitions using the bash script ALIGN_HN.sh.
% 
% Requires COMPUTE_RMSE.m and COMPUTE_XSIM.m from the FANSI-toolbox, which is
% attached as a submodule.
%
% CHANGELOG:
%
% 2024-11-26 (MTC). Original version, with documentation.

clearvars;
close all;

% Utilities functions contained within this repository
addpath(fullfile('.','utils'));
addpath(genpath(fullfile('.','FANSI-toolbox')));

%% Set up Data Set

% Data directory
dir_data = '/media/cherukara/DATA/HN_Repeatability_BIDS/';      % HN Repeatability Study

% Output data directory
dir_save = '/home/cherukara/Documents/Coding/MTC_QSM/Analysis/HNRepeatability_Data/';

% Subject numbers
subs = 1:10;

% Session numbers
sessions = 1:6;

% Method name
str_unwr = 'SEGUE';
str_mask = 'nfe';
str_bkgr = 'VSHARP';
str_susc = 'autoNDI';

% Put the method name together
str_full = strcat('_unwrapped-',str_unwr,'_mask-',str_mask,'_bfr-',str_bkgr,'_susc-',str_susc,'_');
% str_full = strcat('_unwrapped-',str_unwr,'_bfr-',str_bkgr,'_susc-',str_susc,'_');

% Registered method name
str_meth = 'VSHARP';

% Long ROI names
names_roi = {'Thalamus','Caudate Nucleus','Putamen','Globus Pallidus'};

% Lengths
n_subs = length(subs);
n_sess = length(sessions);
n_rois = length(names_roi); 

% Pre-allocate storage arrays
res_rmse = zeros(n_subs,n_sess,n_rois+1);
res_xsim = zeros(n_subs,n_sess,n_rois+1);



%% Loop Through and Load in the Data

for ss = subs

    % Subject name
    if ss < 10
        subname = strcat('sub-0',num2str(ss));
    else
        subname = strcat('sub-',num2str(ss));
    end

    % Load the session 1 brain mask
    arr_mask = double(niftiread(strcat(dir_data,'rawdata/',subname,'/ses-01/anat/',...
                                subname,'_ses-01_desc-brain_mask')));
    vec_mask = arr_mask(:);

    % Load the session 1 ROIs mask
    arr_rois = double(niftiread(strcat(dir_data,'rawdata/',subname,'/ses-01/anat/',...
                                subname,'_ses-01_desc-neckrois_mask')));

    % Load the session 1 QSM data
    arr_orig = niftiread(strcat(dir_data,'derivatives/qsm/',subname,'/ses-01/qsm/',...
                                subname,'_ses-01',str_full,'Chimap'));
    vec_orig = arr_orig(:);


    % Reference to whole-brain average
    arr_orig = arr_orig - mean(vec_orig(vec_mask==1),'omitnan');

    % Test GP value
    vec_gp = arr_orig(arr_rois == 13);
    if mean(vec_gp) < 0
        disp(['Contrast inverted! (',subname,'_ses-01)']);
    end

    for ww = sessions(2:end)

        % Session and scan name
        sesname = strcat('ses-0',num2str(ww));
        dir_raw = strcat(dir_data,'rawdata/',subname,'/',sesname,'/anat/');
        dir_qsm = strcat(dir_data,'derivatives/qsm/',subname,'/',sesname,'/qsm/');
        scanname = strcat(subname,'_',sesname);

        % Load the data from this session
        arr_susc = niftiread(strcat(dir_qsm,scanname,'_mask-brain_reg-ses01_method-',str_meth,'_Chimap'));
        vec_susc = arr_susc(:);

        % Reference
        arr_susc = arr_susc - mean(vec_susc(vec_mask==1),'omitnan');

        % Test GP value
        vec_gp = arr_susc(arr_rois == 13);
        if mean(vec_gp) < 0
            disp(['Contrast inverted! (',scanname,')']);
        end

        % Loop over ROIs and calculate metrics for each one
        for rr = 1:n_rois

            % Pull out ROI mask
            mask_roi = arr_rois == (rr + 9);
            
            % Calculate RMSE 
            res_rmse(ss,ww,rr) = 0.01*compute_rmse(vec_orig(mask_roi(:)),vec_susc(mask_roi(:)));

            % Calculate XSIM
            res_xsim(ss,ww,rr) = compute_xsim(arr_orig,arr_susc,mask_roi);

        end % for rr = 1:n_rois

        % Calculate RMSE and XSIM for the whole brain
        res_rmse(ss,ww,end) = 0.01*compute_rmse(vec_orig(vec_mask==1), vec_susc(vec_mask==1));
        res_xsim(ss,ww,end) = compute_xsim(arr_orig, arr_susc, arr_mask);

    end % for ww = sessions

    fprintf('\tCompleted %s\n',subname);

end % for ss = subs


%% Save the Data
save(strcat(dir_save,'Voxel_Repeatability_',str_meth,'_data.mat'),...
    'res_rmse','res_xsim','str_full','subs','sessions','names_roi');
