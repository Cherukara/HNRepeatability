function [Output,Slices] = ZeroFilling(Input,Dimensions,FirstVoxel)
%DESCRIPTION: Output=ZeroFilling(Input,Dimensions)
%             0 fills Input image up to Dimensions
%
%INPUTS:
%   Input(double matrix) - input image
%   Dimensions(double vector) - output image dimensions
%   FirstVoxel(double vector) - Coordinates of the first pixel in the new,
%                               zero padded image
%
%OUTPUTS:
%   Output(double matrix) - output, 0 filled image
%
%DEPENDENCIES:
%   None
%
%AUTHOR:
%   Anita Karsa, University College London, 2015

% Check Dimensions
if length(size(Input))~=length(Dimensions)
    errordlg('Number of dimensions of Input and length of vector Dimensions has to be equal.');
end

Output = zeros(Dimensions);

% Determine slices
Slices = cell(1,length(Dimensions));
for j = 1:length(Dimensions)
    switch nargin
        case 2
            Slices{j} = floor(Dimensions(j)/2-size(Input,j)/2+1):floor(Dimensions(j)/2+size(Input,j)/2);
        case 3
            Slices{j} = FirstVoxel(j):(FirstVoxel(j)+size(Input,j)-1);
    end
end

% Fill Output with Input
switch length(Dimensions)
    case 1
        Output(Slices{1}) = Input;
    case 2
        Output(Slices{1},Slices{2}) = Input;
    case 3
        Output(Slices{1},Slices{2},Slices{3}) = Input;
    otherwise
        wrndlg('Number of dimensions can be 1,2 or 3');
end



