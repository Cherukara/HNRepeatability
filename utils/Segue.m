function Unwrapped = Segue(Inputs)
%DESCRIPTION: Unwrapped = Segue(Inputs)
%             Speedy rEgion-Growing algorithm for Unwrapping Estimated
%             phase
%             Inputs can be single-echo (3D) or multi-echo (4D)
%
%INPUTS:
%   Inputs.Phase(double array) - Wrapped phase image (3D or 4D) in any units 
%   Inputs.Magnitude(double array) - Magnitude image (3D or 4D), optional 
%                                    (please provide either a magnitude 
%                                    image or a tissue mask)
%   Inputs.Mask(double array) - Mask (3D), optional 
%                               (please provide either a magnitude 
%                               image or a tissue mask)
%   Inputs.Twopi(double) - Voxel value difference corresponding to 2*pi 
%                         (default = 2*pi)
%   Inputs.Preq(double) - Ratio (between 0 and 1) of the entire volume 
%                         unwrapped using a limit of 30% of the border 
%                         (default = 0.7, reduce for faster, increase for 
%                         more accurate unwrapping)
%
%OUTPUTS:
%   Unwrapped(double matrix) - Unwrapped phase image in the original units
%
%DEPENDENCIES:
%   ak_Partitioning_3D_efficient.m
%   ak_SingleStep_3D_efficient.m
%   ak_Neighbourhood_3D.m
%   ak_BorderSizes.m
%
%AUTHOR:
%   Anita Karsa, University College London, 2017

%% Check parameters

Unwrapped = [];

% Phase
if ~isfield(Inputs,'Phase')|~ismember(length(size(Inputs.Phase)),[3,4])|min(size(Inputs.Phase))==1|~or(isa(Inputs.Phase,'double'),isa(Inputs.Phase,'single'))
    warndlg('Please provide a 3D or 4D phase volume (double)!','!! Warning !!');
    return
end

% Mask or Magnitude
if ~isfield(Inputs,'Magnitude')&~isfield(Inputs,'Mask')
    warndlg('Please provide either a tissue mask or the magnitude image!','!! Warning !!');
    return
end

% Masking
switch isfield(Inputs,'Mask')
    case 1
        if sum(size(Inputs.Mask)-size(Inputs.Phase(:,:,:,1)))~=0|~or(isa(Inputs.Mask,'double'),isa(Inputs.Mask,'single'))
            warndlg('Tissue mask has to be double and the same size as the phase image(s)!','!! Warning !!');
            return
        end
    case 0
        if sum(size(Inputs.Magnitude)-size(Inputs.Phase))~=0|~or(isa(Inputs.Magnitude,'double'),isa(Inputs.Magnitude,'single'))
            warndlg('Magnitude image has to be double and the same size as the phase image!','!! Warning !!');
            return
        else %Masking of the magnitude
            switch length(size(Inputs.Magnitude))
                case 3
                    [ftshist, binpos] = hist(Inputs.Magnitude(:),100);
                    
                    diffs = diff(ftshist);
                    diffs = abs(diffs);
                    diffs = diffs./diffs(1);
                    
                    thbin = find(diffs<0.001);
                    thbin = thbin(1);
                    Threshold = binpos(thbin);
                    
                    Inputs.Mask = zeros(size(Inputs.Magnitude));
                    Inputs.Mask(Inputs.Magnitude>=Threshold) = 1;
                case 4
                    Size = size(Inputs.Magnitude);
                    M = reshape(Inputs.Magnitude,[prod(Size(1:3)) Size(4)]);
                    TE = (0:(Size(4)-1))';
                    W = sqrt((sum(M.^2,2))./(sum(M.^2,2).*sum((M.^2)*(TE.^2),2)-sum(M.^2*TE,2).^2));
                    W = reshape(W,Size(1:3));
                    
                    W = abs(W);
                    W(isnan(W)) = 0;
                    W(isinf(W)) = 0;
                    Nstd = 1./W;
                    Nstd(isnan(Nstd)) = 0;
                    Nstd(isinf(Nstd)) = 0;
                    Inputs.Mask = zeros(size(Nstd));
                    Inputs.Mask(Nstd>mean(Nstd(:))*1.1) = 1;
                    
                    clear M W Nstd 
            end
        end
end

clear Inputs.Magnitude

% Twopi
if ~isfield(Inputs,'Twopi')
    Inputs.Twopi = 2*pi;
elseif ~isa(Inputs.Twopi,'double')
    warndlg('Twopi has to be double!','!! Warning !!');
    return
end

% Preq
if ~isfield(Inputs,'Preq')
    Inputs.Preq = 0.7;
elseif ~isa(Inputs.Preq,'double')|Inputs.Preq<0|Inputs.Preq>1
    warndlg('Preq has to be a double between 0 and 1!','!! Warning !!');
    return
end

%% Only unwrap largest connected region
ConnectedRegions = bwlabeln(Inputs.Mask,6);
[C,b] = unique(sort(ConnectedRegions(:)));
b(end+1) = numel(Inputs.Mask);
b = diff(b);
masks = [C,b];
clear C b
masks = masks(masks(:,1)~=0,:);
LabMax = masks(masks(:,2)==max(masks(:,2)),1);

Inputs.Mask(ConnectedRegions~=LabMax) = 0;
clear ConnectedRegions masks LabMax

%% Segue (the main bit)
Unwrapped = zeros(size(Inputs.Phase));

for k = 1:size(Inputs.Phase,4)
    disp(['Echo ',num2str(k)]);
    
    Inputs.Phase(:,:,:,k) = Inputs.Phase(:,:,:,k).*Inputs.Mask;
    
    %% Partitioning
    Labels = ak_Partitioning_3D_efficient(Inputs.Phase(:,:,:,k),Inputs.Mask,Inputs.Twopi);
    
    %% Unwrapping
    [Unwrapped(:,:,:,k),~,Newmask(:,:,:,k)] = ak_SingleStep_3D_efficient(Inputs.Phase(:,:,:,k),Labels,Inputs.Twopi,Inputs.Preq);

end

MaskCommon = prod(Newmask,4);

%% Centering and thresholding the images
for k = 1:size(Inputs.Phase,4)
    Unwrapped_Single = Unwrapped(:,:,:,k);
    Unwrapped_Single = (Unwrapped_Single - Inputs.Twopi*round(mean(Unwrapped_Single(MaskCommon==1))/Inputs.Twopi)).*Newmask(:,:,:,k);
    
    Mean = mean(Unwrapped_Single(Newmask(:,:,:,k)==1));
    SD = std(Unwrapped_Single(Newmask(:,:,:,k)==1));
    Min = Mean - 10*SD;
    Max = Mean + 10*SD;
    Outliers = find(Unwrapped_Single>Max);
    Unwrapped_Single(Outliers) = Unwrapped_Single(Outliers) - Inputs.Twopi*round((Unwrapped_Single(Outliers)-Max)/Inputs.Twopi);
    Outliers = find(Unwrapped_Single<Min);
    Unwrapped_Single(Outliers) = Unwrapped_Single(Outliers) + Inputs.Twopi*round((Min-Unwrapped_Single(Outliers))/Inputs.Twopi);
    
    Unwrapped(:,:,:,k) = Unwrapped_Single;
end



        