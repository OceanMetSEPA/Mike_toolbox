function varargout=plot(meshStruct,varargin)
% Basic plotter for MIKE
%
% INPUT:
% struct with xMesh,yMesh,zMesh coordinates

options=struct;
options.width=900;
options.close=0;
options.data2Plot=[];
%options.maxColVal=[];
options.colourRange=[];
options.concColourMap=flipud(autumn(64));
options.crop=false;
options.dx=0.1;
options.edge='none';
options.boundary='k';
options.label=false;
options.figureHandle=[];
options.status=true;
options.wspz=[];
options.colorbar=true;
options.textMarkers=[];
options.fontsize=8;
options.rotation=0;
options.title='';
options.alpha=0.3;
options.bathy=true;
options.interpreter='tex';
options.zMesh=[];
options.return='figure';
options.ndp=2; % for colour bar
options=checkArguments(options,varargin);


if ~isempty(options.figureHandle) && options.status
    options.status=false;
end

% Sort figure size:
w=options.width;
dx=max(meshStruct.xMesh)-min(meshStruct.xMesh);
dy=max(meshStruct.yMesh)-min(meshStruct.yMesh);
aspectRatio=dy/dx;
figureHeight=aspectRatio*w;
%fprintf('Figure size = %.0f by %.0f (aspect ratio = %f)\n',w,figureHeight,aspectRatio)

% Generate figure:
if isempty(options.figureHandle)
    mainFigureHandle=prepareFigure('close',options.close,'width',w,'height',figureHeight);
else
    mainFigureHandle=options.figureHandle;
end
title(options.title,'fontsize',16,'interpreter',options.interpreter);

% Plot boundary?
if ~isempty(options.boundary)
    try
        plot(meshStruct.xMeshBoundary,meshStruct.yMeshBoundary,'color',options.boundary,'linewidth',3) % domain boundary
    catch
        %        fprintf('Mesh boundary not found\n')
    end
end
% Plot bathymetry?
if options.bathy
    try
        % Colour map for bathy:
        bathyColourMap=gray(64); % B&W so it doesn't interfere with coloured concentrations
        % Get RGB values to represent depth of each cell:
        bathyColours=getColourMatrix(meshStruct.zMesh,bathyColourMap);
        % Now plot bathy triangles:
        zMesh=options.zMesh;
        if isempty(zMesh)
            zMesh=meshStruct.zMesh*0;
        end
        zMesh(isnan(zMesh))=0;
        trisurf(meshStruct.meshIndices,meshStruct.xMesh,meshStruct.yMesh,zMesh,'FaceVertexCData',bathyColours,'EdgeColor',options.edge,'FaceAlpha',0.7);
    catch
        warning('Can''t plot bathy!')
        options.bathy=false;
    end
end

% Plot concs?
data2Plot=options.data2Plot;
if ~isempty(data2Plot)
    if ischar(data2Plot) % Assume it's field of struct
        try
            fn=fieldnames(meshStruct);
            fni=char(closestStringMatch(fn,data2Plot));
            data2Plot=meshStruct.(fni);
            title(fni,'fontsize',14)
        catch
            error('Problem extracting data from struct')
        end
    end
    data2Plot=full(data2Plot);
    data2Plot=data2Plot(:);
    colRange=options.colourRange;
    if isempty(colRange)
        colRange=[min(data2Plot),max(data2Plot)];
    end
    concColourMap=options.concColourMap;
    concColours=getColourMatrix(data2Plot,concColourMap,colRange);
    concAlpha=double(data2Plot>min(colRange));
    z=meshStruct.zMesh*0; % This has nans! Which aren't plotted by trisurf function
    z=zeros(size(z)); % Use this instead
    dataHandle=trisurf(meshStruct.meshIndices,meshStruct.xMesh,meshStruct.yMesh,z,'FaceVertexCData',concColours,'EdgeColor','none','FaceVertexAlphaData',concAlpha,'FaceAlpha','flat');
    if options.crop % Crop axis to non-zero data
        k=data2Plot>0;
        if ~isfield(meshStruct,'x') % Rectangular meshes generated in python don't have x (cell centres) - should maybe fix this sometime...
            meshStruct.x=mean(meshStruct.xMesh(meshStruct.meshIndices),2);
            meshStruct.y=mean(meshStruct.yMesh(meshStruct.meshIndices),2);
        end
        ax=boundaryRectangle(meshStruct.x(k),meshStruct.y(k),'axis',1,'dx',options.dx);
        axis(ax);
    end
    if options.colorbar
        cbTicks=prettyIntervals(colRange);
        drawColourBar(concColourMap,cbTicks);
    end
end

% Plot WSPZ - struct array with Polyshape field
wspz=options.wspz;
if ~isempty(wspz)
    Nwspz=length(wspz);
    pzint=arrayfun(@(x)intersect(wspz(x).Polyshape,meshStruct.boundaryPolyshape),1:Nwspz);
    k=arrayfun(@(x)isempty(x.Vertices),pzint);
    wspz(k)=[];
    Nwspz=length(wspz);
    for i=1:Nwspz
        s=wspz(i);
        plot(s.Polyshape,'facecolor','b','facealpha',options.alpha);
    end
end

% Plot text markers - struct with fields x,y,text (e.g. farm locations)
if ~isempty(options.textMarkers)
    farmLocationStruct=options.textMarkers;
    xFarm=farmLocationStruct.x;
    yFarm=farmLocationStruct.y;
    farmNames=farmLocationStruct.text;
    text(xFarm,yFarm,farmNames,'fontsize',options.fontsize,'rotation',options.rotation)
end
view(2)
grid off
zoom on
% Prepare output:
if nargout==1
    if strcmp(options.return,'figure')
        varargout{1}=mainFigureHandle;
    elseif strcmp(options.return,'data')
        varargout{1}=dataHandle;
    end
end
% Display cell number - probably don't want this (slow and messy) but was
% useful for testing...
if options.label
    k=1:length(meshStruct.x);
    str=arrayfun(@num2str,k,'unif',0);
    text(meshStruct.x(k),meshStruct.y(k),str)
end

set(mainFigureHandle,'Visible','on')

if options.status
    statusbarHandle=statusbar(mainFigureHandle,''); % Need figure to be visible prior to calling statusbar function
    newFont=java.awt.Font('Consolas', java.awt.Font.PLAIN, 12); % Uni-spaced font
    statusbarHandle.getComponent(0).setFont(newFont);
    set(statusbarHandle,'Text','Status bar...')
    set (mainFigureHandle, 'WindowButtonMotionFcn', @mouseMove);
    assignin('base','sb',statusbarHandle)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOUSE STUFF (Move, click)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mouse move
    function mouseMove(~,~)
        % global xMousePos
        % global yMousePos
        % global meshStruct
        % global statusbarHandle
        % global wspz
        % global data2Plot
        % update mouse position
        C=get(gca,'CurrentPoint');
        xMousePos=C(1,1);
        yMousePos=C(1,2);
        meshIndex=Mike.meshIndex(xMousePos,yMousePos,meshStruct);
        try
            try
                depth=meshStruct.z(meshIndex);
            catch
                tri=meshStruct.meshIndices(meshIndex,:);
                depth=mean(meshStruct.zMesh(tri));
            end

            depthString=sprintf('(depth=%.1fm)',depth);
        catch % err
            %    disp(err)
            depthString='';
        end
        statusbarText=sprintf('%.3f, %.3f %s: Cell %d',xMousePos,yMousePos,depthString,meshIndex);
        if ~isempty(wspz)
            id=whichWSPZ(xMousePos,yMousePos,wspz);
            if id
                str=sprintf('Protection Zone : %s (%s)',wspz(id).WSPZ_NAME,tdisp(wspz(id).WSPZ_ID));
            else
                str='Outwith protection zone';
            end
            statusbarText=sprintf('%s           %s',statusbarText,str);
        end
        if ~isempty(data2Plot) && ~isnan(meshIndex)
            statusbarText=sprintf('%s : data val = %f',statusbarText,data2Plot(meshIndex));
        end
        set(statusbarHandle,'Text',statusbarText);
    end

    function id=whichWSPZ(x,y,wspz)
        % Find which wild salmon protection zone (if any) contains coordinates x,y
        %
        % INPUTS:
        % x - longitude
        % y - latitude
        % wspz - struct array loaded from sealiceZones.mat (created by loadSeaLiceProtectionZoneShapefiles.m function)
        %
        % OUTPUT:
        % id - index of wspz struct, or false if coordinates outside
        %
        % EXAMPLE: find water body of every cell in model domain
        % id=whichWSPZ(meshInfo.x,meshInfo.y,wspz); % takes ~ 4 minutes
        %
        Nx=length(x);
        if Nx>1
            id=arrayfun(@(i)whichWSPZ(x(i),y(i),wspz),1:Nx)';
            return
        end

        NZones=length(wspz);
        for zoneIndex=1:NZones
            s=wspz(zoneIndex);
            %    id=inpolygon(x,y,s.Longitude,s.Latitude);
            id=isinterior(s.Polyshape,x,y);
            if id
                id=zoneIndex;
                return
            else
                id=0;
            end
        end

    end

end