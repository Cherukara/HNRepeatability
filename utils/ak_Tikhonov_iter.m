function SusceptibilityMap = ak_Tikhonov_iter(FrequencyMap, Weights, Parameters)
%DESCRIPTION: SusceptibilityMap=ak_Tikhonov_iter(FrequencyMap, Weights, Parameters)
%             Calculates susceptibility map from input image based on with
%             Tikhonov regularisation
%
%INPUTS:
%   FrequencyMap(double matrix) - input image in ppm !!! Needs to be axial (use sag2ax.m or cor2ax.m if necessary)
%   Weights(double matrix) - Weights.*Mask
%   Parameters(struct): Parameters.Alpha(double) - regularisation parameter
%                       Parameters.Resolution(double vector) - image resolution vector (dx,dy,dz) in mm
%                       Parameters.Orientation(double vector) - 3-element vector for defining the z axis 
%                                                               OR
%                                                               6-element vector in the DICOM file (Image orientation or (0020,0037))
%                       Parameters.MatrixSize(double vector) - matrix size before Fourier transformation (0 filled if necessary)
%                       Parameters.SchweserCorr(string) - 'Yes' or 'No'
%                       Parameters.Order - 0 or 1
%
%OUTPUTS:
%   SusceptibilityMap(double matrix) - output image in ppm
%
%DEPENDENCIES:
%   ZeroFilling.m
%   ak_Dk.m
%   ak_ConjGrad.m
%
%AUTHOR: 
%   Anita Karsa, University College London, 2016

%% Define mask and weights
Mask = Weights;
Mask(Mask~=0) = 1;

Weights = 1./Weights;
Weights(isnan(Weights)) = 0;
Weights(isinf(Weights)) = 0;
Weights = Weights/mean(Weights(Mask==1));

%% Define k-space kernel
Parameters.MatrixSize = size(FrequencyMap);
Dk = ak_Dk(Parameters);

%% Define b = DW^2 f
b = real(fftshift(ifftn(ifftshift(Dk.*fftshift(fftn(ifftshift(Weights.^2.*FrequencyMap)))))));
b = b(:);

%% Define system matrix A = DW^2D + alpha*M^2

A=@(xx)(SystemMatrix(Dk,Mask,Weights,Parameters.Alpha,xx));

%% Solve
x = ak_ConjGrad(A,b,Parameters.Threshold,'No');

%% PSF correction (using the direct kernel)

Corr = ifftn(ifftshift(Dk.^2./(Dk.^2 + Parameters.Alpha)));
x = x/Corr(1);

%% Remove background fields
SusceptibilityMap = reshape(x,size(Dk)).*Mask;

end

function b = SystemMatrix(Dk,Mask,Weights,Alpha,xx)

xx = reshape(xx,size(Dk));
b = real(fftshift(ifftn(ifftshift(Dk.*fftshift(fftn(ifftshift(xx)))))));
b = real(fftshift(ifftn(ifftshift(Dk.*fftshift(fftn(ifftshift(Weights.^2.*b)))))));
b = b + Alpha*Mask.*xx;
b = b(:);

end

