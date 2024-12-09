function BorderSizes = ak_BorderSizes(Labels,LinearIndices)
%DESCRIPTION: BorderSizes = ak_StageOne_3D(Labels)
%             Counts voxels on the borders for all regions 
%
%INPUTS:
%   Labels(integer matrix) - Labels
%
%OUTPUTS:
%   BorderSizes(integer matrix) - Lab x 2 containing [labels, border sizes]
%
%DEPENDENCIES:
%   None
%
%AUTHOR: 
%   Anita Karsa, University College London, 2017

% Neighbouring labels of all voxels
Neighbour_Labels = ak_Neighbourhood_3D(LinearIndices,Labels,1);
Neighbour_Labels = Neighbour_Labels(:);
% Labels of all voxels
Labels_Lin = repmat(Labels(LinearIndices),[1 6]);
Labels_Lin = Labels_Lin(:);

RegionBorders = [Labels_Lin, Neighbour_Labels];

clear Labels_Lin Neighbour_Labels;

% Borders are where the labels of the neighbouring pixels are different
RegionBorders = RegionBorders((RegionBorders(:,1)-RegionBorders(:,2))~=0,:);

% We are not interested in Labels == 0
RegionBorders = RegionBorders(RegionBorders(:,1)~=0,:);
RegionBorders = RegionBorders(RegionBorders(:,2)~=0,:);

% We are not interested in the neighbouring labels (we only want to count them)
RegionBorders = RegionBorders(:,1);

% Sort labels to be able to count them
RegionBorders = sort(RegionBorders);

% Identify the type of borders
[C,b] = unique(RegionBorders);
% Count the voxels on each border
b(end+1) = length(RegionBorders)+1;
b = diff(b);

% Labels and border sizes
BorderSizes = [C, b];