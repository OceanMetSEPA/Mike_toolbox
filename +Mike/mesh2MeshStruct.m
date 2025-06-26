function meshStruct=mesh2MeshStruct(f,varargin)
% Get mesh properties / cell areas for coordinates in volumeStruct
%
% This function involves converting grid vertices from lat/lon to
% easting/northing which takes a bit of time and only really needs to be
% done one.
%
% INPUT:
% meshFileName
%
% Optional inputs:
% matfileName [] - filename of matfile to load from / save to
% loadFromMatfile [true] - load from matfile if it exists (faster)
% save [true] - save to matfile
% overwrite [false] - only overwrite existing matfile if this is true
%
% OUTPUT:
% struct containing fields:
%     {'xMesh'        }
%     {'yMesh'        }
%     {'zMesh'        }
%     {'meshIndices'          }
%     {'xVertices'         }
%     {'yVertices'         }
%     {'zVertices'         }
%     {'xMeshBoundary'}
%     {'yMeshBoundary'}
%     {'easting'      }
%     {'northing'     }
%     {'cellArea'     }
%     {'boundaryPolyshape'}
%     {'meshPolyshape'};
%
switch class(f)
    case 'char'
        if ~isfile(f)
            error('Char input should be filename')
        end
        p=fileparts(f);

        volumeStruct=Mike.loadMesh(f,'bound',0);
    case 'struct'
        volumeStruct=f;
        p=pwd; % for now - could add this as an option?
    otherwise
        error('invalid input')
end
matfileName=fullfile(p,'meshInfo.mat');

% Mesh coordinates:
xMesh=volumeStruct.xMesh;
yMesh=volumeStruct.yMesh;
zMesh=volumeStruct.zMesh;
% Triangulation:
meshIndices=volumeStruct.meshIndices;
% Vertex coordinates
xVertices=xMesh(meshIndices);
yVertices=yMesh(meshIndices);
zVertices=zMesh(meshIndices);

volumeStruct.x=mean(xVertices,2);
volumeStruct.y=mean(yVertices,2);

% Get eastings/northings of grid (so we can calculate area of grid cells)
%[easting,northing]=OS.catCoordinates(trix,triy); % Pre 20220516 way - inefficient cos we're doing lots of duplicate vertices!
% !!! mesh might already be in eastings/northings!
if all(yMesh<90)
    [easting,northing]=OS.catCoordinates(xMesh,yMesh);
else
    easting=xMesh;
    northing=yMesh;
end
% Boundary - looks nice in plot :-)
[~,xMeshBoundary,yMeshBoundary]=meshBoundary(meshIndices,xMesh,yMesh);

% Get size info about mesh triangles:
trixEasting=easting(meshIndices);
triyNorthing=northing(meshIndices);
pp=polygonProperties(trixEasting,triyNorthing);
%cellArea=pp.area; % Use this rather than volume/dz from volumeStruct as it is constant in time

% Prepare output struct
meshStruct=struct;
meshStruct.x=volumeStruct.x;
meshStruct.y=volumeStruct.y;
meshStruct.xMesh=xMesh;
meshStruct.yMesh=yMesh;
meshStruct.zMesh=zMesh;
meshStruct.meshIndices=meshIndices;
meshStruct.xVertices=xVertices;
meshStruct.yVertices=yVertices;
meshStruct.zVertices=zVertices;
meshStruct.xMeshBoundary=xMeshBoundary;
meshStruct.yMeshBoundary=yMeshBoundary;
meshStruct.easting=easting;
meshStruct.northing=northing;

meshStruct.cellArea=pp.area;%
%meshStruct.cellProperties=pp;
meshStruct.boundaryPolyshape=polyshape(meshStruct.xMeshBoundary,meshStruct.yMeshBoundary,'simplify',0);

% Also generate polyshape of individiual elements. This takes a while (few
% minutes) but is useful for checking intersections for generateing tracks
Nt=length(meshStruct.xVertices);
nans=nan(Nt,1);
xt=[meshStruct.xVertices,nans]';
yt=[meshStruct.yVertices,nans]';
xt=xt(:);
yt=yt(:);
ps=polyshape(xt,yt,'simplify',0);
meshStruct.meshPolyshape=ps;

var2matfile(meshStruct,matfileName)

end

