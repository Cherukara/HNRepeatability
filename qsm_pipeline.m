% MTC_QSM_PIPELINE.m - A QSM processing pipeline for all applications
%
% This script is based on my previous "MTC_xx_PIPELINE.m" scripts, but designed
% to be unified across all applications, and take a properly labelled and
% categorized BIDS dataset
%
% It is designed to run for one data set at a time, rather than to loop across
% many.
%
%
%       Copyright (C) University College London, 2024
%
%
% Created by MT Cherukara, May 2024
%
% CHANGELOG:
%
% 2024-06-10 (MTC). Changed the way the masking is implemented so that processes
%       a series of non-exclusive masking steps (brain mask, noise-based, MFG,
%       and connected components) and then logically ANDs the mask
%
% 2024-10-08 (MTC). Refactored a little bit so that everything is in the correct
%       units throughout.
%
% 2024-10-24 (MTC). Added functionality to calculate the brain mask (using BET)
%       within this script, rather than relying on an external BET mask.
%
% 2024-10-30 (MTC). Added a bit more versatility to the SWITCH-CASE sections

clearvars;

% mtc_setup;

%% Select the Data Set

% Data directory
dir_data = '/media/cherukara/DATA/HERD_Study_BIDS/';                 % HERD Study
% dir_data = '/media/cherukara/DATA/HN_Repeatability_BIDS/';      % HN Repeatability Study
% dir_data = '/media/cherukara/DATA/DECOMPOSE/';
% dir_data = 'D:\Matthew\7_Academia\Data\Perspectum_SS\';  % Perspectum Study, COMPUTWO
% dir_data = '/media/cherukara/DATA/Perspectum_SS_BIDS/';         % Perspectum Study
% dir_data = '/media/cherukara/DATA/SourceSeparation_Study_BIDS/';    % Source-Sep Study

% Subject number
for sub = 4

% Session number - set this to 0 if there is no "session" level in the BIDS
% hierarchy
for ses = 1

% Do we have the echoes as separate niftis (as per BIDS format)
is_echoes = 1;

% Key Parameters (change these by hand if necessary)
Params.Orientation = [0, 0, 1];
Params.Vsize = 22;      % V-SHARP kernel size
Params.Alpha = 0.06;        % IterTik regularization parameter


%% Select Processing Options

% Choose algorithm for each step
method_fitting = 'load';
method_unwrap = 'SEGUE';
method_bgfr = 'PDF';
method_dipole = 'autoNDI';

% Optional Processing Steps
is_MPPCA = 1;
is_MSMV = 0;

% Masking steps
is_mask_brain = 0;
is_mask_noise = 1;
is_mask_mfg = 0;
is_mask_conn = 1;
is_mask_erode = 1;

% Choose which outputs to save
save_mask = 1;
save_unwrap = 1;
save_localfield = 1;

% Start the timer
tStart = tic;



%% Define Scan Name

% Determine subject name
if sub < 10
    subname = strcat('sub-0',num2str(sub));
else
    subname = strcat('sub-',num2str(sub));
end

% Set the folders and define SCANNAME
if ses == 0
    dir_raw = strcat(dir_data,'rawdata/',subname,'/anat/');
    dir_qsm = strcat(dir_data,'derivatives/qsm/',subname,'/qsm/');
    scanname = subname;
else
    sesname = strcat('ses-0',num2str(ses));
    dir_raw = strcat(dir_data,'rawdata/',subname,'/',sesname,'/anat/');
    dir_qsm = strcat(dir_data,'derivatives/qsm/',subname,'/',sesname,'/qsm/');
    scanname = strcat(subname,'_',sesname);
end

% Define method names
if strcmp(method_unwrap,'load')
    method_unwrap_name = 'SEGUE';
else
    method_unwrap_name = method_unwrap;
end

if strcmp(method_bgfr,'load')
    method_bgfr_name = 'VSHARP';
else
    method_bgfr_name = method_bgfr;
end

if is_MSMV == 1
    method_bgfr_name = strcat(method_bgfr_name,'mSMV');
end

% if is_mask_brain == 1
%     method_unwrap_name = strcat(method_unwrap_name,'brain');
% end


%% Read basic information
% Attempt to load a JSON file from the rawdata directory and read from it some
% key information, if this doesn't work, we'll just assume some default values

try 
    if is_echoes == 1
        txt_json = fileread(strcat(dir_raw,scanname,'_part-mag_echo-1_GRE.json'));
    else
        txt_json = fileread(strcat(dir_raw,scanname,'_part-mag_GRE.json'));
    end
    data_json = jsondecode(txt_json);

    
catch
    disp('Loading JSON did not work');
    data_json = struct([]);
end

% Read some basic information from the JSON file
Params = mtc_parse_json(data_json,Params);


% Read in information from the NIFTI
if is_echoes == 1
    inf_mag = niftiinfo(strcat(dir_raw,scanname,'_part-mag_echo-1_GRE'));
else
    inf_mag = niftiinfo(strcat(dir_raw,scanname,'_part-mag_GRE'));
end

Params.Resolution = inf_mag.PixelDimensions(1:3);
Params.MatrixSize = inf_mag.ImageSize(1:3);
Params.NEchoes = 5;


%% Load the Data


if is_echoes == 1

    % Pre-allocate MAG and PHA arrays
    arr_mag = zeros([Params.MatrixSize,Params.NEchoes]);
    arr_pha = zeros([Params.MatrixSize,Params.NEchoes]);

    % Loop through and load magnitude and phase
    for ee = 1:Params.NEchoes
    
        % Read in Echo times
        txt_json = fileread(strcat(dir_raw,scanname,'_part-mag_echo-',num2str(ee),'_GRE.json'));
        data_json = jsondecode(txt_json);
        Params.TEs(ee) = data_json.EchoTime;


    
        % Load magnitude and phase
        arr_mag(:,:,:,ee) = double(niftiread(strcat(dir_raw,scanname,'_part-mag_echo-',num2str(ee),'_GRE')));
        arr_pha(:,:,:,ee) = double(niftiread(strcat(dir_raw,scanname,'_part-phase_echo-',num2str(ee),'_GRE')));
    
    end

else

    % Load the data
    arr_mag = double(niftiread(strcat(dir_raw,scanname,'_part-mag_GRE')));
    arr_pha = double(niftiread(strcat(dir_raw,scanname,'_part-phase_GRE')));

    % Number of echo times
    Params.NEchoes = size(arr_mag,4);

end

% Manually specify the echo times
%       Eventually, this should be read in from a .json in the BIDS folder
Params.TEs = (1:Params.NEchoes).*4.92e-3;

% Remove NaNs
arr_mag(isnan(arr_mag)) = 0;
arr_pha(isnan(arr_pha)) = 0;

% Scale the phase between 0 and 2*pi
arr_pha = arr_pha - min(arr_pha(:));
arr_pha = 2*pi*arr_pha./max(arr_pha(:));

% Calculate delta TE and the Phase scaling
if ~isfield(Params,'CF')
    Params.CF = 123.239;        % Carrier Frequency (in MHz)
end

Params.dTE = Params.TEs(2) - Params.TEs(1);
Params.rad2Hz = 2*pi*Params.dTE;
Params.rad2ppm = 2*pi.*Params.TEs(1).*Params.CF;


%% 1. Non-Linear Fitting of Complex Data

switch lower(method_fitting)

    case {'medi','fitppm'}

        fprintf('Calculating Non-Linear Fit... \n');
        tLocal = tic;

        % Make complex data
        arr_comp = arr_mag.*exp(-1i*arr_pha);

        % MP-PCA denoising
        if is_MPPCA == 1
            % See if denoising has already been done, if so, load that data
            try
                arr_mag = niftiread(strcat(dir_raw,scanname,'_part-mag_denoised-CV_GRE'));
                arr_pha = niftiread(strcat(dir_raw,scanname,'_part-phase_denoised-CV_GRE'));

                arr_comp = arr_mag.*exp(-1i*arr_pha);

            % Otherwise, calculate it here and save it!
            catch

                fprintf('\tRunning MP-PCA denoising... \n');
                arr_comp = denoiseCV(arr_comp,[2,2,2]);
    
                % Save out the denoised data
                arr_mag = abs(arr_comp);
                arr_pha = angle(arr_comp);
    
                % Update info struct so that it saves better
                inf_mag.ImageSize = size(arr_mag);
                inf_mag.PixelDimensions = [Params.Resolution, 1];
                inf_mag.Datatype = 'double';
                inf_mag.BitsPerPixel = 64;
    
                % Save denoised multi-echo magnitude and phase
                niftiwrite(arr_mag, strcat(dir_raw,scanname,'_part-mag_denoised-CV_GRE'), inf_mag);
                niftiwrite(arr_pha, strcat(dir_raw,scanname,'_part-phase_denoised-CV_GRE'), inf_mag);
            end
        end

        % Non-linear fitting using MEDI toolbox function
        [arr_fieldfit, arr_noise, ~, arr_phizero] = Fit_ppm_complex_TE(arr_comp,Params.TEs);

        % Absolutize the noise
        arr_noise = abs(arr_noise);

        % Calculate an R2* map
        arr_R2s = arlo(Params.TEs,arr_mag);

        % Create an INFO struct to save the Fit and Noise data
        inf_noise = inf_mag;
        inf_noise.ImageSize = size(arr_noise);
        inf_noise.PixelDimensions = Params.Resolution;
        inf_noise.Datatype = 'double';
        inf_noise.BitsPerPixel = 64;

        % Save Fit and Noise
        niftiwrite(arr_fieldfit, strcat(dir_qsm,scanname,'_acq-MEDI_fit'), inf_noise);
        niftiwrite(arr_noise, strcat(dir_qsm,scanname,'_acq-MEDI_noise'), inf_noise);
        niftiwrite(arr_R2s, strcat(dir_qsm,scanname,'_acq-ARLO_R2smap'), inf_noise);

        fprintf('\tCompleted Non-Linear Fitting \n');
        tPart = toc(tLocal);
        fprintf('\t\tTime taken: %.4f \n', tPart);

    case 'load'

        % Load the noise
        arr_noise = niftiread(strcat(dir_qsm,scanname,'_acq-MEDI_noise'));

        % Load the INFO struct from the noise data
        inf_noise = niftiinfo(strcat(dir_qsm,scanname,'_acq-MEDI_noise'));

        % Load the fitted field
        arr_fieldfit = niftiread(strcat(dir_qsm,scanname,'_acq-MEDI_fit'));

        % Load the R2* map
        arr_R2s = niftiread(strcat(dir_qsm,scanname,'_acq-ARLO_R2smap'));

end

%% 2. Masking

% For masking, the steps in the process are additive, rather than exclusive, so
% instead of a switch-case construction we just use a set of consecutive IF
% statements and multiply the masks together

% Start with a universal mask
arr_mask = ones(size(arr_noise));

% Mask description string
str_mask = '_desc-';

% Brain masking
if is_mask_brain == 1

    try
        % Try to load a brain mask
        arr_tempmask = niftiread(strcat(dir_raw,scanname,'_desc-brain_mask'));

    catch
        % If not, calculate it manually by running BET
        fprintf('\tCalculating BET mask... \n');
        arr_tempmask = BET(squeeze(arr_mag(:,:,:,1)),Params.MatrixSize,Params.Resolution);
        
        % Save out the newly discovered BET mask
        inf_mask = inf_mag;
        inf_mask.ImageSize = size(arr_mask);
        inf_mask.PixelDimensions = Params.Resolution;
        inf_mask.Datatype = 'int16';
        inf_mask.BitsPerPixel = 16;
        niftiwrite( int16(arr_tempmask), strcat(dir_raw,scanname,'_desc-brain_mask'), inf_mask);

    end


    % Apply the brain mask to the main mask
    arr_mask = arr_mask .* double(arr_tempmask);

    % Update description string
    str_mask = strcat(str_mask,'b');

end 

% Noise-based masking
if is_mask_noise == 1

    % Create noise-based mask
    Weights = 1./arr_noise;
    Weights(isnan(Weights)) = 0;
    Weights(isinf(Weights)) = 0;
    arr_tempmask = Weights > mean(Weights(:))*1.2;

    % Apply the current mask to the main mask
    arr_mask = arr_mask .* arr_tempmask;

    % Update description string
    str_mask = strcat(str_mask,'n');

end

% Mask based on magnitude of field gradients
if is_mask_mfg == 1

    % Gradient-based threshold mask
    arr_tempmask = GradientBasedThreshold(arr_fieldfit, arr_mask, 3) == 1;

    % Apply the current mask to the main mask
    arr_mask = arr_mask .* arr_tempmask;

    % Update description string
    str_mask = strcat(str_mask,'g');

end

% Mask out unconnected components of the original mask
if is_mask_conn == 1

    CC = bwconncomp(arr_mask,18);
    Careas = regionprops(CC, 'Area');
    LL = labelmatrix(CC);
    arr_tempmask = ismember(LL, find([Careas.Area] >= 0.5*max([Careas.Area])));

    % Apply the current mask to the main mask
    arr_mask = arr_mask .* arr_tempmask;

    % Update description string
    str_mask = strcat(str_mask,'c');

end

if is_mask_erode == 1

    % Erode the mask using a radius 1 disk
    arr_tempmask = imerode(arr_mask,strel('disk',1));

    % Apply the current mask to the main mask
    arr_mask = arr_mask .* arr_tempmask;

    % Update mask description string
    str_mask = strcat(str_mask,'e');

end


% Save the mask
if save_mask == 1
    inf_mask = inf_mag;
    inf_mask.ImageSize = size(arr_mask);
    inf_mask.PixelDimensions = Params.Resolution;
    inf_mask.Datatype = 'int16';
    inf_mask.BitsPerPixel = 16;
    niftiwrite( int16(arr_mask), strcat(dir_qsm,scanname,str_mask,'_mask'), inf_mask);
end


%% 3. Phase Unwrapping

fprintf('Unwrapping Phase... \n');
tLocal = tic;

switch lower(method_unwrap)

    % Laplacian unwrapping using AK's function which calls the MEDI toolbox
    case {'laplace','laplacian','lpu'}

        Params.MatrixSize = [512 512 512];
        arr_field = ak_Laplacian_MEDI(arr_fieldfit, Params);
        Params.MatrixSize = size(arr_fieldfit);


    % SEGUE unwrapping using AK's function
    case 'segue'
        SParams.Phase = arr_fieldfit;
        SParams.Mask = double(arr_mask);
        arr_field = Segue(SParams);


    % Load a pre-calculated SEGUE unwrapped fieldmap
    case 'load'

        % Don't save it again
        save_unwrap = 0;

        % Load NIFTI
        arr_field = niftiread(strcat(dir_qsm,scanname,'_part-phase',...
                                     '_unwrapped-',method_unwrap_name,'_field'));

end

% Save the unwrapped field map (in Radians)
if save_unwrap == 1

    % Basic info structure
    inf_unwrap = inf_mag;
    inf_unwrap.ImageSize = size(arr_field);
    inf_unwrap.PixelDimensions = Params.Resolution;
    inf_unwrap.Datatype = 'double';
    inf_unwrap.BitsPerPixel = 64;

    % Save
    niftiwrite(double(arr_field), ...
               strcat(dir_qsm,scanname,'_part-phase_unwrapped-',method_unwrap,...
               '_field'), inf_unwrap);
end

tPart = toc(tLocal);
fprintf('\t\tTime taken: %.4f \n', tPart);


%% 4. Background Field Removal

fprintf('Removing Background Fields... \n');
tLocal = tic;

switch lower(method_bgfr)

    case 'pdf'

        % PDF takes arr_field in RADIANS and returns a arr_fieldloc in RADIANS

        % Perform PDF background field removal using MEDI function
        arr_fieldloc = PDF(arr_field, arr_noise, arr_mask, ...
                        Params.MatrixSize, Params.Resolution, Params.Orientation, 6e-2 );

        % Apply the mask to the local field
        arr_fieldloc = arr_fieldloc .* arr_mask;

        % SMV Prefilter on
        Params.SMVPrefilter = 1;

    case {'vsharp','v-sharp'}

        % VSHARP takes input field in RADIANS and returns local field in RADIANS

        % Perform V-SHARP (using STI Suite)
        [arr_fieldloc, arr_mask] = V_SHARP(arr_field, arr_mask,...
                                   'voxelsize',Params.Resolution,...
                                   'smvsize',Params.Vsize);

        % SMV Prefilter off
        Params.SMVPrefilter = 0;

        % Save a newly eroded V-SHARP mask
        if save_mask == 1
            inf_mask = inf_mag;
            inf_mask.ImageSize = size(arr_mask);
            inf_mask.PixelDimensions = Params.Resolution;
            inf_mask.Datatype = 'int16';
            inf_mask.BitsPerPixel = 16;
            niftiwrite( int16(arr_mask), strcat(dir_qsm,scanname,str_mask,'VS_mask'), inf_mask);
        end


    case {'iharp','iharperella','harperella','harp'}

        % iHARP takes input field in RADIANS and returns local field in RADIANS

        % Calculate IHARPERELLA (using STI Suite)
        arr_fieldloc = iHARPERELLA(arr_field, arr_mask, 'voxelsize', Params.Resolution);

        % SMV Prefilter on
        Params.SMVPrefilter = 1;

    case 'nlpdf'

        % nlPDF takes input field in RADIANS and returns local field in RADIANS

        % Create a dipole kernel
        PParams.K = dipole_kernel_fansi(Params.MatrixSize, Params.Resolution, 1);

        % Input fieldmap in RADIANS
        PParams.input = arr_field;

        % Other inputs
        PParams.mask = arr_mask;
        PParams.weight = arr_mask.*max(arr_mag,[],4);
        PParams.spatial_res = Params.Resolution;

        % Calculate non-linear PDF
        pdf_out = npdfCG(PParams);

        % Store output
        arr_fieldloc = double(pdf_out.local).*arr_mask;

    case 'tfi'

        % If we are going to use TFI for the dipole inversion, then we shouldn't
        % do anything at this stage
        arr_fieldloc = arr_field;
        save_localfield = 0;
        is_MSMV = 0;
        method_dipole = 'TFI';

    case 'load'

        % Don't save it again
        save_localfield = 0;

        % Load previously calculated V-SHARP data
        arr_fieldloc = niftiread(strcat(dir_qsm,scanname,...
                              '_unwrapped-',method_unwrap_name,...
                              '_bfr-',method_bgfr_name,'_localfield'));

end

% mSMV filtering
if is_MSMV == 1
    arr_fieldloc = msmv(arr_fieldloc, arr_mask, arr_R2s, Params.Resolution, 5, 5, 10, ...
                     Params.B0, Params.SMVPrefilter);
end

% Save the local field
if save_localfield == 1

    % Basic info structure
    inf_local = inf_mag;
    inf_local.ImageSize = size(arr_fieldloc);
    inf_local.PixelDimensions = Params.Resolution;
    inf_local.Datatype = 'double';
    inf_local.BitsPerPixel = 64;

    % Save
    niftiwrite(arr_fieldloc, strcat(dir_qsm,scanname,...
               '_unwrapped-',method_unwrap_name,...
               '_bfr-',method_bgfr_name,'_localfield'),...
               inf_local);

end

fprintf('\tCompleted Backgrond Field Removal \n');
tPart = toc(tLocal);
fprintf('\t\tTime taken: %.4f \n', tPart);


%% 5. Dipole Inversion

fprintf('Inverting Dipole... \n');
tLocal = tic;

% At this point, arr_fieldloc should be in RADIANS

% Generic FANSI parameters
FParams.voxelSize = Params.Resolution;
FParams.maxOuterIter = 1000;
FParams.tolUpdate = 0.02;

% Generic FANSI data - what units should the fieldmap be in??
FParams.input = arr_fieldloc;
FParams.weight = max(arr_mag,[],4)./max(arr_mag(:));
FParams.mask = arr_mask;
FParams.K = dipole_kernel_fansi(Params.MatrixSize, Params.Resolution, 0);

switch lower(method_dipole)

    case {'dirtik','direct'}

        % Requires input in PPM and produces output in PPM
        arr_fieldloc = arr_fieldloc ./ Params.rad2ppm;

        % Direct Tikhonov (using Anita Karsa's code)
        Params.Alpha = 0.11;
        arr_susc = ak_Tikhonov(arr_fieldloc.*arr_mask, arr_mask, Params);

    case {'itertik','tikhonov','tik'}

        % Requires input in PPM and produces output in PPM
        arr_fieldloc = arr_fieldloc ./ Params.rad2ppm;

        % Iterative Tikhonov (using Anita Karsa's code)
        Params.Threshold = 6e-2;
        arr_susc = ak_Tikhonov_iter(arr_fieldloc.*arr_mask, arr_noise.*arr_mask, Params);

    case 'tfi'

        % NOT SURE ABOUT THE UNITS!!

        % Parameters
        Params.Field = arr_fieldloc;
        Params.Magnitude = max(arr_mag,[],4);
        Params.Mask = arr_mask;
        Params.Noise = arr_noise;
        Params.Lambda = 1e-6;
        Params.R2s = arr_R2s;

        % Do TFI (using my code)
        arr_susc = mtc_tfi(Params);

    case {'fansi','nltv','nl-tv'}

         % Fansi Options
         FParams.alpha1 = 1e-5;
         FParams.mu1 = 100.*FParams.alpha1;
         FParams.mu2 = 1.0;

         % Run FANSI
         struct_fansi = nlTV(FParams);
         arr_susc = struct_fansi.x.*arr_mask;

    case {'whfansi','wh-fansi','wh-nltv','wh-nl-tv'}

        % Fansi Options
        FParams.alpha1 = 2e-4;
        FParams.beta = 1e4;
        FParams.mu1 = 100.*FParams.alpha1;
        FParams.mu2 = 1.0;

        % Run FANSI
        struct_fansi = WH_nlTV(FParams);
        arr_susc = struct_fansi.x;

    case {'starqsm','star-qsm','qsm-star','qsmstar'}

        % StarQSM takes an input fieldmap in RADIANS and returns a
        % susceptibility distribution in PPB 

        % Run Star QSM from STI Suite
        arr_susc = QSM_star(arr_fieldloc, arr_mask,...
            'TE',Params.TEs(1),'B0',Params.B0,...
            'H',Params.Orientation,'padsize',[12,12,12],...
            'voxelsize',Params.Resolution);

        % Rescale into PPM from PPB
        arr_susc = arr_susc./1000;

    case 'ilqsr'

        % iLSQR takes an input fieldmap in RADIANS and returns a susceptibility
        % distribution in PPB

        % Run iLSQR from STI Suite
        arr_susc = QSM_iLSQR(arr_fieldloc, arr_mask,...
            'TE',Params.TEs(1),'B0',Params.B0,...
            'H',Params.Orientation,'padsize',[12,12,12],...
            'voxelsize',Params.Resolution);

        % Rescale into PPM from PPB
        arr_susc = arr_susc./1000;

    case {'l1qsm','l1','wl1tv'}

        % Fansi Options
        FParams.alpha1 = 1e-3;
        FParams.lambda = 10;
        FParams.mu1 = 125.*FParams.alpha1;
        FParams.mu2 = 1.0;

        % FANSI inputs
        FParams.weight = FParams.lambda.*FParams.weight;

        % Run FANSI
        struct_fansi = wL1TV(FParams);
        arr_susc = struct_fansi.x;

    case {'ndi','autondi'}

        % autoNDI definitely takes input fieldmap in RADIANS and returns a
        % susceptibility map in RADIANS

        % FANSI Options
        FParams.input = arr_fieldloc;
        FParams.alpha = 1e-2; % regularization weight
        FParams.tau = 1; % gradient descent rate

        % Run FANSI
        struct_fansi = ndi_auto(FParams);

        % Rescale the result from RADIANS into PPM 
        arr_susc = struct_fansi.x.*arr_mask./Params.rad2ppm;

    case {'medi'}

        % NOT SURE ABOUT THE UNITS

        % MEDI parameters
        temp_file_MEDI = 'tmp_RDF.mat';
        
        % Everything has to be named exactly right
        iFreq = arr_field.*Params.PhaScale;      % Unwrapped phase, in radians
        RDF = arr_fieldloc.*Params.PhaScale;        % Fieldmap, in radians
        N_std = arr_noise;      % Noise standard deviation
        iMag = max(arr_mag,[],4);       % Magnitude image 
        Mask = arr_mask;       % Mask of ROI
        matrix_size = Params.MatrixSize;
        voxel_size = Params.Resolution;
        delta_TE = Params.dTE;
        CF = Params.CF;
        B0_dir = Params.Orientation;


        % Save a .mat file in the current folder
        save(temp_file_MEDI,'iFreq','RDF','N_std','iMag','Mask','matrix_size','voxel_size','delta_TE','CF','B0_dir');

        % Run MEDI
        tic;
        arr_susc = MEDI_L1_2('filename',temp_file_MEDI,'lambda',1000,'merit',...
                             'smv',5,'data_weighting',1,'gradient_weighting',1);
        toc;

        % % This is the MEDI_L1 call from the SEPIA toolbox:
        % chi = MEDI_L1('filename',tmp_filename,'lambda',lambda,'data_weighting',wData,'gradient_weighting',wGrad,...
        %       'merit',isMerit,'smv',radius,'zeropad',pad,'lambda_CSF',lam_CSF,'percentage',percentage);

        % Delete the temporary file
        delete(temp_file_MEDI);


    case 'qsmnet'

        % Chi Sep Tool directory
        dir_chisep = '/home/cherukara/Documents/Coding/Toolboxes/chi-separation/old/Chi_Sep_Tool';

        % Convert the input fieldmap from RADIANS to HZ
        arr_fieldhz = arr_fieldloc./Params.rad2Hz;

        % For large shapes, we have to split the array into two slabs:
        n_slabs = 3;
        width_slab = Params.MatrixSize(3)./n_slabs;

        % Pre-allocate result array
        arr_susc = zeros(Params.MatrixSize);

        for ss = 1:n_slabs

            i_slab = (ss-1).*width_slab + 1;
            e_slab = ss.*width_slab;

            % Define a slab and zero-pad
            slab_fieldhz = padarray(arr_fieldhz(:,:,i_slab:e_slab),[0,0,2]);
            slab_mask = padarray(arr_mask(:,:,i_slab:e_slab),[0,0,2]);

            % Do QSMNet
            slab_susc = QSMnet_general(dir_chisep, slab_fieldhz, slab_mask, slab_mask, ...
                                      Params.Orientation, Params.CF, ...
                                      Params.Resolution, Params.MatrixSize );

            % Put the slabs together (undoing zero-padding)
            arr_susc(:,:,i_slab:e_slab) = slab_susc(:,:,3:(end-2));

        end

end

% Create info structure
inf_susc = inf_mag;
inf_susc.ImageSize = size(arr_susc);
inf_susc.PixelDimensions = Params.Resolution;
inf_susc.Datatype = 'double';
inf_susc.BitsPerPixel = 64;

% Save the susceptibility map
niftiwrite(double(arr_susc), strcat(dir_qsm,scanname,...
           '_unwrapped-',method_unwrap_name,...
           '_bfr-',method_bgfr_name,...
           '_susc-',method_dipole,'_Chimap'), inf_susc);

% Timings
fprintf('\tCompleted Dipole Inversion \n');
tPart = toc(tLocal);
fprintf('\t\tTime taken: %.4f \n', tPart);



%% Wrap Up

% Timings
fprintf('Completed %s\n',scanname);
tEnd = toc(tStart);
fprintf('Total time elapsed: %.4f \n', tEnd);

end
end
