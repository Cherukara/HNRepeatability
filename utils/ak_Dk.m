function [Dk,K2] = ak_Dk(Parameters)
%DESCRIPTION: Dk = ak_Dk(Parameters)
%             Calculates k-space kernel for forward calculation
%
%INPUTS:
%   Parameters(struct): Parameters.Resolution(double vector) - image resolution vector (dx,dy,dz) in mm
%                       Parameters.Orientation(double vector) - 3-element vector for defining the z axis 
%                                                               OR
%                                                               6-element vector in the DICOM file (Image orientation or (0020,0037))
%                       Parameters.MatrixSize(double vector) - matrix size before Fourier transformation (0 filled if necessary)
%
%OUTPUTS:
%   Dk(double matrix) - Kernel
%
%DEPENDENCIES:
%   ZeroFilling.m
%
%AUTHOR: 
%   Anita Karsa, University College London, 2015

% Sort parameters

%Resolution
Resolution = Parameters.Resolution;
%MatrixSize
if isfield(Parameters,'MatrixSize')
    MatrixSize = Parameters.MatrixSize;
else
    warndlg('Define matrix size');
end
%Image orientation (z_prime = real B0 direction in our coordinate system)
if isfield(Parameters,'Orientation')
    switch length(Parameters.Orientation)
        case 6
            z_prime = cross(Parameters.Orientation(1:3),Parameters.Orientation(4:6));
            z_prime = z_prime/norm(z_prime);
        case 3
            z_prime = Parameters.Orientation;
            z_prime = z_prime/norm(z_prime);
        otherwise
            z_prime = [0 0 1];
    end
else
    z_prime = [0 0 1];
end

% Create dipole kernel

dkx = 1/Resolution(1)/MatrixSize(1);
dky = 1/Resolution(2)/MatrixSize(2);
dkz = 1/Resolution(3)/MatrixSize(3);

Center=ceil(MatrixSize/2);

Range = cell(3,1);

for Index = 1:3
    if mod(MatrixSize(Index),2)==0
        Range{Index} = -Center(Index):Center(Index)-1;
    else
        Range{Index} = -Center(Index)+1:Center(Index)-1;
    end
end

[Y,X,Z] = meshgrid(Range{2}*dky,...
                   Range{1}*dkx,...
                   Range{3}*dkz);

K2 = X.^2+Y.^2+Z.^2;
               
Dk = 1/3 - (z_prime(1)*X + z_prime(2)*Y + z_prime(3)*Z).^2./K2;
Dk(isnan(Dk)) = 0;

