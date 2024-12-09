function Neighbours = ak_Neighbourhood_3D(LinearIndices,Matrix,Spacing)
%DESCRIPTION: Neighbours = ak_Neighbourhood_3D(LinearIndices,Matrix,Spacing)
%             Finds the "Spacing" order neighbours in 3D of each voxel in 
%             LinearIndeces and returns their value
%
%INPUTS:
%   LinearIndices(double vector) - Indices of voxels whose neighbours we
%                                  are interested in
%   Matrix(double array) - Array with the values
%   Spacing(integer) - Order of neighbours
%
%OUTPUTS:
%   Neighbours(double array, Nx6) - Values at each neighbour
%
%DEPENDENCIES:
%   None
%
%AUTHOR: 
%   Anita Karsa, University College London, 2016

Size = size(Matrix);

[X,Y,Z] = ind2sub(Size, LinearIndices);

Neighbours = zeros(length(LinearIndices),6);

% X-
X_ = X;
X_(X_>Spacing) = X_(X_>Spacing)-Spacing;
Neighbours(:,1) = Matrix(sub2ind(Size,X_,Y,Z));

% X+
X_ = X;
X_(X_<(Size(1)-Spacing+1)) = X_(X_<(Size(1)-Spacing+1))+Spacing;
Neighbours(:,2) = Matrix(sub2ind(Size,X_,Y,Z));

% Y-
Y_ = Y;
Y_(Y_>Spacing) = Y_(Y_>Spacing)-Spacing;
Neighbours(:,3) = Matrix(sub2ind(Size,X,Y_,Z));

% Y+
Y_ = Y;
Y_(Y_<(Size(2)-Spacing+1)) = Y_(Y_<(Size(2)-Spacing+1))+Spacing;
Neighbours(:,4) = Matrix(sub2ind(Size,X,Y_,Z));

% Z-
Z_ = Z;
Z_(Z_>Spacing) = Z_(Z_>Spacing)-Spacing;
Neighbours(:,5) = Matrix(sub2ind(Size,X,Y,Z_));
 
% Z+
Z_ = Z;
Z_(Z_<(Size(3)-Spacing+1)) = Z_(Z_<(Size(3)-Spacing+1))+Spacing;
Neighbours(:,6) = Matrix(sub2ind(Size,X,Y,Z_));