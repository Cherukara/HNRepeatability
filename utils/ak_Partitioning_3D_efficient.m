function Labels = ak_Partitioning_3D_efficient(WrappedImage, TissueMask, Twopi)
%DESCRIPTION: Labels = ak_Partitioning_3D_efficient(WrappedImage, TissueMask)
%             Performs partitioning similar to Jenkinson, MRM, 2003
%             2pi interval is diveded into Nintervals smaller intervals
%             Phase values in each interval are partitioned after erosion
%             Eroded voxels are assigned to regions
%
%INPUTS:
%   WrappedImage(double matrix) - Wrapped image in any units
%   TissueMask(double matrix) - binary tissue mask
%   Nintervals(positive integer) - Number of smaller intervals 2pi is split
%                                  into
%
%OUTPUTS:
%   Labels(double matrix) - Labels
%
%DEPENDENCIES:
%   ak_Neighbourhood_3D.m
%
%AUTHOR: 
%   Anita Karsa, University College London, 2017

% Mask image
WrappedImage = WrappedImage.*TissueMask;

% Determine the smaller intervals
Min = min(WrappedImage(:));
Max = max(WrappedImage(:));
Step = Twopi/6;
Borders = Min:Step:(Max+0.01);

% Set 'background' to a value outside the borders so it wouldn't be included
% in the partitioning
WrappedImage(TissueMask==0) = Max+1;

% Initialise Labels array
Labels = zeros(size(WrappedImage));

Size = size(Labels);

% h = waitbar(0,'Partitioning in 3D...');

for j = 1:(length(Borders)-1) % loop across each smaller interval
    
    % Find voxels within range
    mask = zeros(size(WrappedImage));
    mask(WrappedImage>=Borders(j) & WrappedImage<Borders(j+1)) = 1;

    % Erode mask where it is thin in at least 2 directions -> mask2
    mask2 = mask;
    LinearIndices = bwperim(mask,6);
    LinearIndices = find(LinearIndices==1);
    FirstOrder = ak_Neighbourhood_3D(LinearIndices,1-mask,1);
    SecondOrder = ak_Neighbourhood_3D(LinearIndices,1-mask,2);
    ThirdOrder = ak_Neighbourhood_3D(LinearIndices,1-mask,3);
    OppositeSide = [FirstOrder(:,2),FirstOrder(:,1),FirstOrder(:,4),FirstOrder(:,3),...
        FirstOrder(:,6),FirstOrder(:,5)];
    % IsThere = +1 where: 0 1 0 or 0 1 1 0 or 0 1 1 1 0
    IsThere = FirstOrder.*OppositeSide + SecondOrder.*OppositeSide + ThirdOrder.*OppositeSide;
    clear FirstOrder SecondOrder ThirdOrder FourthOrder OppositeSide 
    IsThere = [IsThere(:,1)+IsThere(:,2), IsThere(:,3)+IsThere(:,4),...
               IsThere(:,5)+IsThere(:,6)];
    IsThere(IsThere~=0) = 1;
    IsThere = IsThere(:,1)+IsThere(:,2)+IsThere(:,3);
    mask2(LinearIndices(IsThere>=2)) = 0;
    clear IsThere    
    Edges = mask - mask2;
    
    % Assign labels in the eroded mask (mask2)
    Labels_ = bwlabeln(mask2,6);
    
    Edges2 = zeros(size(mask2)+2);
    Edges2(2:(end-1),2:(end-1),2:(end-1)) = mask2;
    Edges2 = bwperim(1-Edges2,6);
    Edges2 = Edges2(2:(end-1),2:(end-1),2:(end-1))+mask2;
    
    Neighbours = zeros(3,3,3);
    Neighbours([5 11 13 14 15 17 23]) = 1;

    StoppingCriterium = 1;
    while StoppingCriterium~=0
      
        Indices = find(Edges.*Edges2==1);
        % Edges: 1. Look at neighbouring voxels
        Neighbour_Labels = ak_Neighbourhood_3D(Indices,Labels_,1);
        
        Neighbour_Labels = sort(Neighbour_Labels,2);
        Labels_(Indices) = Neighbour_Labels(:,end);
        
        [X,Y,Z] = ind2sub(size(mask2),Indices);
        Edges2(Indices) = 1;
        X_ = X; X_(X_>1) = X_(X_>1)-1; Edges2(sub2ind(Size,X_,Y,Z)) = 1;
        X_ = X; X_(X_<Size(1)) = X_(X_<Size(1))+1; Edges2(sub2ind(Size,X_,Y,Z)) = 1;
        Y_ = Y; Y_(Y_>1) = Y_(Y_>1)-1; Edges2(sub2ind(Size,X,Y_,Z)) = 1;
        Y_ = Y; Y_(Y_<Size(2)) = Y_(Y_<Size(2))+1; Edges2(sub2ind(Size,X,Y_,Z)) = 1;
        Z_ = Z; Z_(Z_>1) = Z_(Z_>1)-1; Edges2(sub2ind(Size,X,Y,Z_)) = 1;
        Z_ = Z; Z_(Z_<Size(3)) = Z_(Z_<Size(3))+1; Edges2(sub2ind(Size,X,Y,Z_)) = 1;
               
        %Edges2 = Edges2 + convn(Edges.*Edges2,Neighbours,'same');
        %Edges2(Edges2>1) = 1;
        
        Edges(Indices) = 0;
        
        StoppingCriterium = sum(Neighbour_Labels(:,end));
    end
    
    % Update Labels
    Edges(Labels_~=0) = 0;
    mask(Edges==1) = 0;
    Labels = Labels + mask.*(Labels_ + max(Labels(:)));
    
    % Edges: 2. Assign labels to the rest (small, disconnected regions of a 
    % few voxels, where every neighbour was region 0)   
    Labels_ = bwlabeln(Edges,6);
    Labels = Labels + Edges.*(Labels_ + max(Labels(:)));
    
    % waitbar(j/(length(Borders)-1),h);
end

% close(h);




