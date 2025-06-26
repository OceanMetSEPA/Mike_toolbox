function [op, cells2Keep] = cropMesh(meshStruct,cells2Keep)
% Crop struct containing mesh information
%
% INPUTS:
%  meshStruct - struct containing xMesh, meshIndices etc.
% Either:
% 'ax' - keep points within axis
% 'polyshape' - keep points within this polyshape
% 'indices' - keep these indices
%
% NB this function is slow if polyshape is a complicated shape and there
% are lots of points to query
%
%
if islogical(cells2Keep)
    if isequal(size(cells2Keep,1),size(meshStruct.meshIndices,1))
        % default - we've specified logical values for each cell
    elseif isequal(size(cells2Keep,1),size(meshStruct.xMesh,1))
        k=find(cells2Keep);
        cells2Keep=any(ismember(meshStruct.meshIndices,k),2);
        size(cells2Keep)
    else
        error('Need logical value for each grid cell')
    end
else
    if isa(cells2Keep,'polyshape')
        [xb,yb]=polyshape2polygon(cells2Keep);
    elseif isa(cells2Keep,'double')
        try 
             [xb,yb]=axis2polygon(cells2Keep);
        catch err
            disp(err)
            error('Expected axis for numeric input')
        end
    else
        error('invalid input')
    end
    % Extract mesh data
    meshIndices = meshStruct.meshIndices;
    xv = meshStruct.xMesh(meshIndices);
    yv = meshStruct.yMesh(meshIndices);

    % Identify cells to keep
    cells2Keep = any(inpolygon(xv, yv, xb, yb), 2);
end

if ~any(cells2Keep)
    op = [];
    return;
end

% Extract triangles and renumber them efficiently
trik = meshStruct.meshIndices(cells2Keep, :);
[ind, ~, nind] = unique(trik(:)); % Get unique points and new indices
ntrik = reshape(nind, size(trik)); % Reshape back into triangle format

% Identify points to keep
Nv = numel(meshStruct.xMesh);
points2Keep = false(1, Nv);
points2Keep(ind) = true;

% Filter and update mesh struct
op = structFilter(meshStruct, cells2Keep);
op = structFilter(op, points2Keep);
op.meshIndices = ntrik;
[boundaryIndices,xMeshBoundary,yMeshBoundary]=meshBoundary(op);
op.boundaryIndices=boundaryIndices;
op.xMeshBoundary=xMeshBoundary;
op.yMeshBoundary=yMeshBoundary;
end
