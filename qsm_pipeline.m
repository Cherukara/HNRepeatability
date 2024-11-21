% QSM_PIPELINE.m reconstructs QSMs from BIDS-formatted NIFTI data
%
%
%       Copyright (C) University College London, 2024
%
%
% Created by MT Cherukara, May 2024

clearvars;

% Add toolboxes to path - replace these with your own path
addpath(genpath('/home/cherukara/Documents/Coding/Toolboxes/MEDI_toolbox/'));
addpath(genpath('/home/cherukara/Documents/Coding/Toolboxes/STISuite_V3.0/'));

% Utilities functions contained within this repository
addpath(fullfile('.','utils'));
addpath(genpath(fullfile('.','FANSI-toolbox')));
addpath(genpath(fullfile('.','chi-separation')));


%% Select the Data Set

% Data directory
dir_data = '/media/cherukara/DATA/HN_Repeatability_BIDS/';      % HN Repeatability Study


% Key Parameters (change these by hand if necessary)
Params.Orientation = [0, 0, 1];
Params.Vsize = 22;      % V-SHARP kernel size
Params.Alpha = 0.06;        % IterTik regularization parameter

% Hard-coded parameters (eventually, these should be read from a JSON file in
% the BIDS directory)
Params.B0 = 3;
Params.CF = 123.138;
Params.NEchoes = 4;
Params.TEs = (1:Params.NEchoes).*4.61e-3;


%% Select Processing Options

% Choose algorithm for each step
method_fitting = 'load';
method_unwrap = 'LPU';
method_bgfr = 'PDF';
method_dipole = 'iterTik';

% Optional Processing Steps
is_MPPCA = 1;

% Masking steps
is_mask_noise = 1;
is_mask_mfg = 0;
is_mask_fill = 0;
is_mask_conn = 0;
is_mask_erode = 0;
is_mask_vs = 0;

% Choose which outputs to save
save_mask = 1;
save_unwrap = 1;
save_localfield = 1;


%% Loop Over Subjects and Sessions

% Subject number
for sub = 2:10

% Session number - set this to 0 if there is no "session" level in the BIDS
% hierarchy
for ses = 1:6

% Start the timer
tStart = tic;

%% Define Scan Name

% Determine subject name
if sub < 10
    subname = strcat('sub-0',num2str(sub));
else
    subname = strcat('sub-',num2str(sub));
end

sesname = strcat('ses-0',num2str(ses));
dir_raw = strcat(dir_data,'rawdata/',subname,'/',sesname,'/anat/');
dir_qsm = strcat(dir_data,'derivatives/qsm/',subname,'/',sesname,'/qsm/');
scanname = strcat(subname,'_',sesname);

% Define method names
if strcmp(method_unwrap,'load')
    method_unwrap_name = 'SEGUE';
else
    method_unwrap_name = method_unwrap;
end

if strcmp(method_bgfr,'load')
    method_bgfr_name = 'PDF';
else
    method_bgfr_name = method_bgfr;
end

% Read in information from the NIFTI
inf_mag = niftiinfo(strcat(dir_raw,scanname,'_part-mag_echo-1_GRE'));
Params.Resolution = inf_mag.PixelDimensions(1:3);
Params.MatrixSize = inf_mag.ImageSize(1:3);


%% Load the Data


% Pre-allocate MAG and PHA arrays
arr_mag = zeros([Params.MatrixSize,Params.NEchoes]);
arr_pha = zeros([Params.MatrixSize,Params.NEchoes]);

% Loop through and load magnitude and phase
for ee = 1:Params.NEchoes

    % Do we want to use MP-PCA denoising?
    if is_MPPCA
        try 
            % Load already-denoised Mag and Phase data, if it exists
            arr_mag(:,:,:,ee) = double(niftiread(strcat(dir_raw,scanname,'_part-mag_echo-',num2str(ee),'_denoised-CV_GRE')));
            arr_pha(:,:,:,ee) = double(niftiread(strcat(dir_raw,scanname,'_part-phase_echo-',num2str(ee),'_denoised-CV_GRE')));

            % Don't do the denoising again later
            is_MPPCA = 0;
        catch
            % Load noisy Mag and Phase data
            arr_mag(:,:,:,ee) = double(niftiread(strcat(dir_raw,scanname,'_part-mag_echo-',num2str(ee),'_GRE')));
            arr_pha(:,:,:,ee) = double(niftiread(strcat(dir_raw,scanname,'_part-phase_echo-',num2str(ee),'_GRE')));
        end

    else

        % Load magnitude and phase
        arr_mag(:,:,:,ee) = double(niftiread(strcat(dir_raw,scanname,'_part-mag_echo-',num2str(ee),'_GRE')));
        arr_pha(:,:,:,ee) = double(niftiread(strcat(dir_raw,scanname,'_part-phase_echo-',num2str(ee),'_GRE')));
    end

end

% Remove NaNs
arr_mag(isnan(arr_mag)) = 0;
arr_pha(isnan(arr_pha)) = 0;

% Scale the phase between 0 and 2*pi
arr_pha = arr_pha - min(arr_pha(:));
arr_pha = 2*pi*arr_pha./max(arr_pha(:));

% Calculate echo spacing and scaling parameters
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

        % Non-linear fitting using MEDI toolbox function
        [arr_fieldfit, arr_noise, ~, ~] = Fit_ppm_complex_TE(arr_comp,Params.TEs);

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
str_mask = '';

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

% Fill holes in the mask
if is_mask_fill == 1

    % Imfill the mask
    arr_mask = imfill(arr_mask,'holes');

    % Update description string
    str_mask = strcat(str_mask,'f');

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


if any(strcmp(method_bgfr_name,{'vsharp','v-sharp'}))

    % Update the mask description string to incorporate v-sharp
    str_mask = strcat(str_mask,'v');

    % Don't erode
    is_mask_erode = 0;
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
    niftiwrite( int16(arr_mask), strcat(dir_qsm,scanname,'_desc-',str_mask,'_mask'),...
                inf_mask, 'Compressed',true);
end

% Load the VS mask
if is_mask_vs == 1
    arr_mask = double(niftiread(strcat(dir_qsm,scanname,'_desc-nfv_mask')));
    str_mask = 'nfv';
end


%% 3. Phase Unwrapping

fprintf('Unwrapping Phase... \n');
tLocal = tic;

switch lower(method_unwrap)

    % Laplacian unwrapping using a MEDI toolbox function
    case {'laplace','laplacian','lpu'}
        arr_field = unwrapLaplacian(arr_fieldfit, Params.MatrixSize, Params.Resolution);
        
    % SEGUE unwrapping using AK's function
    case {'segue','exact'}
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
               '_field'), inf_unwrap, 'Compressed',true);
end

tPart = toc(tLocal);
fprintf('\t\tTime taken: %.4f \n', tPart);


%% 4. Background Field Removal

fprintf('Removing Background Fields... \n');
tLocal = tic;

switch lower(method_bgfr)

    case 'pdf'


        % Perform PDF background field removal using MEDI function
        arr_fieldloc = PDF(arr_field, arr_noise, arr_mask, ...
                        Params.MatrixSize, Params.Resolution, Params.Orientation, 6e-2 );

        % Apply the mask to the local field
        arr_fieldloc = arr_fieldloc .* arr_mask;


    case 'lbv'

        % Perform LBV background field removal using MEDI function
        arr_fieldloc = LBV(arr_field, arr_mask, Params.MatrixSize, Params.Resolution);

        % Apply the mask to the local field
        arr_fieldloc = arr_fieldloc .* arr_mask;

    case {'vsharp','v-sharp'}

        % VSHARP takes input field in RADIANS and returns local field in RADIANS

        % Perform V-SHARP (using STI Suite)
        [arr_fieldloc, arr_mask] = V_SHARP(arr_field, arr_mask,...
                                   'voxelsize',Params.Resolution,...
                                   'smvsize',Params.Vsize);


        % Save a newly eroded V-SHARP mask
        if save_mask == 1
            inf_mask = inf_mag;
            inf_mask.ImageSize = size(arr_mask);
            inf_mask.PixelDimensions = Params.Resolution;
            inf_mask.Datatype = 'int16';
            inf_mask.BitsPerPixel = 16;
            niftiwrite( int16(arr_mask), strcat(dir_qsm,scanname,'_desc-',str_mask,'v_mask'),...
                        inf_mask, 'Compressed',true);
        end

    case 'tfi'

        % If we are going to use TFI for the dipole inversion, then we shouldn't
        % do anything at this stage
        arr_fieldloc = arr_field;
        save_localfield = 0;
        method_dipole = 'TFI';

    case 'load'

        % Don't save it again
        save_localfield = 0;

        % Load previously calculated local field data
        arr_fieldloc = niftiread(strcat(dir_qsm,scanname,...
                              '_unwrapped-',method_unwrap_name,...
                              '_bfr-',method_bgfr_name,'_localfield'));

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
               '_mask-',str_mask,...
               '_bfr-',method_bgfr_name,'_localfield'),...
               inf_local, 'Compressed',true);

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
           '_mask-',str_mask,...
           '_bfr-',method_bgfr_name,...
           '_susc-',method_dipole,'_Chimap'), inf_susc,...
           'Compressed',true);

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
