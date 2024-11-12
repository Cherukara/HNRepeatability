function [Unwrapped,NewLabels,Newmask] = ak_SingleStep_3D_efficient(WrappedImage, Labels, Twopi, Preq)
%DESCRIPTION: [Unwrapped,NewLabels,Newmask] = ak_SingleStep_3D_efficient(WrappedImage, Labels, Twopi, Preq)
%             Performs a region-growing unwrapping and merging algorithm
%             which gradually unwraps all neighbours of the largest region
%
%INPUTS:
%   WrappedImage(double matrix) - Wrapped image in any units
%   Labels(integer matrix) - Labels
%   Twopi(double) - Voxel value difference correcponding to 2*pi
%   Preq(double) - Percentage of the entire volume unwrapped using a limit
%                  of 30% of the border
%
%OUTPUTS:
%   Unwrapped(double matrix) - Output image in the original units
%   NewLabels(integer matrix) - New labels
%   NewMask(integer matrix) - Binary tissue mask of unwrapped regions
%
%DEPENDENCIES:
%   ak_Neighbourhood_3D.m
%   ak_BorderSizes.m
%
%AUTHOR:
%   Anita Karsa, University College London, 2017

Plimit = 0.3;

% Initialise neighbourhoods
NeighbourValues = ak_Neighbourhood_3D(1:numel(Labels),WrappedImage,1);
VoxelValues = repmat(WrappedImage(:),[1 6]);
OppositeValues = NeighbourValues(:,[2 1 4 3 6 5]);
Differences = 2*VoxelValues(:)-NeighbourValues(:)-OppositeValues(:);
clear VoxelValues NeighbourValues OppositeValues;

NeighbourLabels = ak_Neighbourhood_3D(1:numel(Labels),Labels,1);
VoxelLabels = repmat(Labels(:),[1 6]);
OppositeLabels = NeighbourLabels(:,[2 1 4 3 6 5]);

% Label1, Label2, Label3, 2*Value2 - Value1 - Value3 (Extrapolation - Actual)
AllNeighbours = [OppositeLabels(:),VoxelLabels(:),NeighbourLabels(:),Differences];
clear OppositeLabels VoxelLabels NeighbourLabels Differences
AllNeighbours = AllNeighbours(AllNeighbours(:,2)~=0,:);
AllNeighbours = AllNeighbours(AllNeighbours(:,3)~=0,:);
AllNeighbours = AllNeighbours((AllNeighbours(:,2)-AllNeighbours(:,3))~=0,:);

% Initialise BorderSizes
BorderSizes = AllNeighbours(:,2);
[C,b] = unique(sort(BorderSizes));
b(end+1) = length(BorderSizes)+1;
BorderSizes = [C diff(b)]; % Labels BorderSize
clear C b
CurrentBorders = BorderSizes;

AllNeighbours = AllNeighbours(AllNeighbours(:,1)~=0,:);

% Initialise Shifts and Sizes
[Shifts,Sizes] = unique(sort(Labels(:)));
Sizes(end+1) = numel(Labels)+1;
Sizes = diff(Sizes);

Sizes = [Shifts Sizes]; % Label, Size
CurrentSizes = Sizes;

NumberOfRegions = length(Shifts);
Shifts = [Shifts zeros(size(Shifts)) zeros(size(Shifts)) Shifts]; % Labels Shift MaskIndicator UpdatedLabel

% Initialise TotalWrapped and TotalVolume
TotalVolume = length(find(Labels~=0));
TotalWrapped = TotalVolume;
TotalWrapped_previous = 0;

NewLabels = Labels;

NumberOfRegionsRest = NumberOfRegions;
% h = waitbar(1/NumberOfRegions,'Regions unwrapped...');

Indicator = 0;

while and(and(TotalWrapped/TotalVolume>0.01,TotalWrapped~=TotalWrapped_previous),Indicator<2)
    
    % Finding the region with the largest border
    if and(TotalWrapped/TotalVolume>(1-Preq),Indicator==0)
        LabMax = CurrentBorders(CurrentBorders(:,2)==max(CurrentBorders(:,2)),1);
        LabMax = LabMax(1);
    elseif and(TotalWrapped/TotalVolume<(1-Preq),Indicator==0)
        LabMax = CurrentSizes(CurrentSizes(:,2)==max(CurrentSizes(CurrentSizes(:,1)~=0,2)),1);
        LabMax = LabMax(1);
        Plimit = 0.1;
        Shifts(:,3) = 0;
        Indicator = 1;
    elseif Indicator==1
        LabMax = CurrentSizes(CurrentSizes(:,2)==max(CurrentSizes(CurrentSizes(:,1)~=0,2)),1);
        LabMax = LabMax(1);
        Plimit = 0;
        Shifts(:,3) = 0;
        Indicator = 2;
    end
    
    if Plimit==0
        TotalWrapped_previous = TotalWrapped;
    end
    
    Shifts(Shifts(:,4)==LabMax,3) = 1;
    
    % Unwrapping and merging neighbouring regions
    Unwr = 0;
    Regions = 0;
    
    while ~or(isempty(Unwr),isempty(Regions))
        
        % Extrapolation of the phase to the neighbours
        Neighbour_Labels = AllNeighbours(and(AllNeighbours(:,1)==LabMax,AllNeighbours(:,2)==LabMax),:);
        Unwr = round(Neighbour_Labels(:,4)/Twopi);
        Neighbour_Labels = Neighbour_Labels(:,3);
        
        if ~isempty(Unwr)
            
            % Selecting regions to unwrap and choosing the corresponding phase shifts
            Regions = ak_SelectRegions(Unwr,Neighbour_Labels,BorderSizes,Plimit); % Votes, Labels, Shift
            
            % Unwrapping of remaining regions
            if ~isempty(Regions)
                
                % Identify different Unwr values
                Unwr = unique(Regions(:,3));
                
                % Apply each phase shift separately
                for k = 1:length(Unwr)
                    
                    % Find Regions with shift Unwr
                    Labs = Regions(Regions(:,3)==Unwr(k),2);
                    
                    % Update all tables
                    Indices = ismember(Shifts(:,4),Labs);
                    Shifts(Indices,2) = Shifts(Indices,2) + Unwr(k);
%                     Region = ismember(Labels,Labs);
%                     WrappedImage(Region) = WrappedImage(Region) + Unwr(k)*Twopi;
%                     Labels(Region) = LabMax;
                    Shifts(Indices,3) = 1;
                    Shifts(Indices,4) = LabMax;
                    
%                     Indices = ismember(Sizes(:,1),Labs);
%                     Sizes(Sizes(:,1)==LabMax,2) = Sizes(Sizes(:,1)==LabMax,2) + sum(Sizes(Indices,2)); 
%                     Sizes = Sizes(~Indices,:); 
                    
                    BorderSizes = BorderSizes(~ismember(BorderSizes(:,1),Labs),:);
                    
                    CurrentBorders = CurrentBorders(~ismember(CurrentBorders(:,1),Labs),:);
                    
                    CurrentSizes = CurrentSizes(~ismember(CurrentSizes(:,1),Labs),:);
                    
                    Indices = ismember(AllNeighbours(:,3),Labs);
                    Done = Indices.*(AllNeighbours(:,1)==LabMax).*(AllNeighbours(:,2)==LabMax);
                    AllNeighbours = AllNeighbours(~Done,:);
                    
                    Indices = ismember(AllNeighbours(:,1),Labs);
                    AllNeighbours(Indices,1) = LabMax;
                    AllNeighbours(Indices,4) = AllNeighbours(Indices,4) - Unwr(k)*Twopi;
                    Indices = ismember(AllNeighbours(:,2),Labs);
                    AllNeighbours(Indices,2) = LabMax;
                    AllNeighbours(Indices,4) = AllNeighbours(Indices,4) + 2*Unwr(k)*Twopi;
                    Indices = ismember(AllNeighbours(:,3),Labs);
                    AllNeighbours(Indices,3) = LabMax;
                    AllNeighbours(Indices,4) = AllNeighbours(Indices,4) - Unwr(k)*Twopi;
                    AllNeighbours = AllNeighbours((AllNeighbours(:,2)-AllNeighbours(:,3))~=0,:);
                    
                    NumberOfRegionsRest = NumberOfRegionsRest - length(Labs);
                    
                end
            end
        end
        
        % waitbar(1-NumberOfRegionsRest/NumberOfRegions,h);
        
    end
    
    % Update borders
    CurrentBorders = CurrentBorders(CurrentBorders(:,1)~=LabMax,:);
    Newmask = zeros(size(Labels));
    Newmask(ismember(Labels,Shifts(Shifts(:,4)==LabMax,1))) = 1;
    NewLabels(Newmask==1) = LabMax;
    CurrentSizes(CurrentSizes(:,1)==LabMax,2) = sum(Newmask(:));
    BorderSizesNew = ak_BorderSizes(NewLabels,find(bwperim(Newmask,6)==1));
    if ~isempty(BorderSizesNew)
        BorderSizes(BorderSizes(:,1)==LabMax,:) = BorderSizesNew;
    end
    
    TotalWrapped = TotalVolume - sum(Sizes(Shifts(:,3)==1,2));
    
end

% close(h);

% Calculate Unwrapped, NewLabels and Newmask

Unwr = unique(Shifts(:,2));
Unwr = Unwr(Unwr~=0);

for k = 1:length(Unwr)
    Labs = Shifts(Shifts(:,2)==Unwr(k),1);
    Region = ismember(Labels,Labs);
    WrappedImage(Region) = WrappedImage(Region) + Unwr(k)*Twopi;
end

% Newmask = zeros(size(WrappedImage));
% Newmask(ismember(Labels,Shifts(Shifts(:,3)==1,1))) = 1;

% NewLabels = Labels;
% NewLabels(Newmask==1) = LabMax;
Unwrapped = WrappedImage.*Newmask;

end

function Regions = ak_SelectRegions(Unwr,Neighbour_Labels,BorderSizes,Plimit)

% Regions = [Labels, Unwr values] for each neighbouring voxel
Regions = [Neighbour_Labels,Unwr];
% Sort rows to gather the same labels + unwr values next to each
% other
Regions = sortrows(Regions);
% Identify neighbouring regions (labels) and candidate phase shifts
% (unwr values)
[C,b] = unique(Regions,'rows');
% Count the voxels for each pair of labels + unwr values
b(end+1) = length(Unwr)+1;
b = diff(b);

% Regions = [Number of votes, Labels, Unwr values]
Regions = [b,C];
% Sort row by Labels
Regions = sortrows(Regions,[2 -1]);
% Find highest votes for each neighbouring region
Diff = diff(Regions(:,2));
Diff = find(Diff~=0);
Diff = [1;Diff+1];
Regions = Regions(Diff,:);

% Count voxels on each border
Neighbour_Labels = sort(Neighbour_Labels);
[C,b] = unique(Neighbour_Labels);
b(end+1) = length(Neighbour_Labels)+1;
b = diff(b);
% Neighbour_Labels = [Border size, Labels, Total border size for neighbour]
Neighbour_Labels = [b,C,BorderSizes(ismember(BorderSizes(:,1),C),2)];

% Check that Regions(:,2) == Neighbour_Labels(:,2), Labels match
if sum(Regions(:,2) - Neighbour_Labels(:,2)) == 0
    % Calculate percentage of votes
    Pagree = Regions(:,1)./Neighbour_Labels(:,1);
    % Calculate percentage of voxels "detected"
    Pborder = Neighbour_Labels(:,1)./Neighbour_Labels(:,3);
else
    disp('Labels do not match!');
    exit;
end

% Keep only the regions where Pagree > f(Pborder)
A = 1/(Plimit-1);
Regions = Regions(Pagree>(A.*(Pborder-1)),:);

end
