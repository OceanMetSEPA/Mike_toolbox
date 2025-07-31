function [gridIndex] = meshIndex(xp, yp,varargin)
% Find the containing quadrilateral or triangle for each point (xp, yp)
%
% INPUT:
% - xp, yp: Query points' coordinates (Nx1 arrays)
% And either:
%    - struct containing xMesh,yMesh,meshIndices
%    - faces,vertices
%    - faces,xMesh,yMesh
%
% where:    
%       faces: (Mx3 or Mx4) list of triangles or quadrilaterals
%       xMesh, yMesh: Mesh node coordinates (Px1 arrays)
%
% OUTPUT:
% - gridIndex: Index of the quadrilateral or triangle containing each (xp, yp) point
%
% Supports **both triangles and quadrilaterals**.
% Uses **tsearchn** for fast point location.

if nargin<3
    help Mike.meshIndex
    return
end

switch length(varargin)
    case 1 % Assume meshStruct input (as returned by Mike.loadMesh for example)
        s=varargin{1};
        try
            fn=fieldnames(s);
            xfield=fn{contains(fn,'xmesh','ig',1)};
            xMesh=s.(xfield);
            yfield=fn{contains(fn,'ymesh','ig',1)};
            yMesh=s.(yfield);
            triField=fn{contains(fn,'triMesh','ig',1)|contains(fn,'MeshIndices','ig',1)|contains(fn,'tri')};
            faces=s.(triField);
            vertices=[xMesh,yMesh];
        catch
            disp(fn)
            error('Can''t find mesh info from argument 3')
        end
    case 2 % faces, vertices
        faces=varargin{1};
        vertices=varargin{2};
        if size(vertices,2)==3
            vertices(:,3)=[];
        end
    case 3 % faces, xMesh, yMesh (original function call)
        faces=varargin{1};
        xMesh=varargin{2};
        yMesh=varargin{3};
        vertices=[xMesh,yMesh];
    otherwise
        error('Incorrect number of arguments')
end

origShape = size(xp);
xp = xp(:);
yp = yp(:);

% Replace nans with Mike's null val
nullVal=Mike.nullVal;
xp(isnan(xp))=nullVal;
yp(isnan(yp))=nullVal;

quadMesh = (size(faces, 2) == 4);  % Check if input is quadrilateral mesh

if quadMesh
    NQuad = size(faces, 1);
    triFaces=[faces(:,1:3);faces(:,[1,3,4])];
else
    triFaces = faces;  % Already triangles
end

%  **Find Containing Triangle**
try
gridIndex = tsearchn(vertices, triFaces, [xp, yp]);  % Returns NaN if outside
catch err
%    fprintf('Point %f, %f\n',xp,yp)
    disp(size(vertices))
    disp(size(triFaces))
    fprintf('Orig shape = %s\n',tdisp(origShape))
    disp(err)
    underline
    gridIndex=nan;
end
%  **Convert Triangle Index Back to Quadrilateral Index (if needed)**
if quadMesh
    gridIndex=mod(gridIndex-1,NQuad)+1;
end

%  **Reshape Output to Match Input**
gridIndex = reshape(gridIndex, origShape);

end
