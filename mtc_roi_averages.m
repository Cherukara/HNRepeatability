% mtc_roi_averages.m
%
% Pull out and store the average values in several ROIs from some Head and Neck
% data

clearvars;
close all;

%% Set up Data Set

% Data directory
dir_data = '/media/cherukara/DATA/HN_Repeatability_BIDS/';      % HN Repeatability Study

% Subject numbers
subs = 1:10;

% Session numbers
sessions = 1:6;

% Long ROI names
names_roi = {'Lymph Node','Lymph Node','Lymph Node','Lymph Node',...
             'Lymph Node','Lymph Node','Parotid Gland','Submandibular Gland'...
             'Thalamus','Caudate Nucleus','Putamen','Globus Pallidus',...
             'Jugular Vein','Trapezius Muscle','Subcutaneous Fat'};

% Method name
str_unwr = 'SEGUE';
str_bkgr = 'PDF';
str_susc = 'StarQSM';

% Should we reference?
is_reference = 1;

% Put the method name together
str_meth = strcat('_unwrapped-',str_unwr,'_bfr-',str_bkgr,'_susc-',str_susc,'_');

% Numbers
n_rois = length(names_roi);
n_subs = length(subs);
n_sess = length(sessions);

% Pre-allocate data arrays
arr_roi_av = zeros(n_subs,n_sess,n_rois);
arr_roi_sd = zeros(n_subs,n_sess,n_rois);


%% Loop Through and Load in Data

for ss = subs

    for ww = sessions

        % Subject name
        if ss < 10
            subname = strcat('sub-0',num2str(ss));
        else
            subname = strcat('sub-',num2str(ss));
        end
    
        % Session and scan name
        sesname = strcat('ses-0',num2str(ww));
        dir_raw = strcat(dir_data,'rawdata/',subname,'/',sesname,'/anat/');
        dir_qsm = strcat(dir_data,'derivatives/qsm/',subname,'/',sesname,'/qsm/');
        scanname = strcat(subname,'_',sesname);
    
        % Load ROI data
        arr_rois = niftiread(strcat(dir_raw,scanname,'_desc-neckrois_mask'));
        vec_rois = arr_rois(:);

        % Load in the data from the QSM-only method
        arr_susc = niftiread(strcat(dir_qsm,scanname,str_meth,'Chimap'));
        vec_susc = arr_susc(:);

        % Referencing?
        if is_reference == 1

            % Load in the reference region mask
            arr_refmask = double(niftiread(strcat(dir_raw,scanname,'_desc-brain_mask')));

            % Calculate reference value
            s_refmean = mean(arr_susc.*arr_refmask,'all','omitnan');

            % Subtract reference value
            arr_susc = arr_susc - s_refmean;

        end % if is_reference == 1

        % Loop over ROIs and store the averages
        for rr = 1:n_rois

            % Pull out the data from the current ROI
            roi_susc = vec_susc(vec_rois == (rr + 1));
    
            % Store the averages
            arr_roi_av(ss,ww,rr) = mean(roi_susc,'omitnan');
            arr_roi_sd(ss,ww,rr) = std(roi_susc,'omitnan');

        end % for rr = 1:n_rois

    end % for ww = sessions

    fprintf('\tCompleted %s\n',subname);

end % for ss = subs


%% Save the Data
dir_save = '/home/cherukara/Documents/Coding/MTC_QSM/Analysis/Analysis_Data/';
save(strcat(dir_save,'Repeatability',str_meth,'ref_data.mat'),...
    'arr_roi_av','arr_roi_sd','str_meth','subs','sessions','names_roi');
