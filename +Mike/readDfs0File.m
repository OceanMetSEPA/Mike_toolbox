function [dfs0Output ] = readDfs0File( dfs0FileName,varargin )
% Extract data from MIKE Dfs0 file
%
% This function extracts time-series data from a MIKE Dfs0 file
% (corresponding to a single point)
%
% INPUT:
% dfs0FileName - file used by MIKE with extension .dfs0
%
% Optional Inputs:
% info (false) - return information about file rather than data (e.g.
%                data parameters, time period, location names)
% parameters [] - return specified parameters
% locations []  - return data for specified location
%
% OUTPUT:
% struct containing metadata about Dfs0 file and data in various formats,
% depending on optional inputs above. Let the number of locations be Nloc
% and the number of parameters by Nparam.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   readDfs0File.m  $
% $Revision:   1.2  $
% $Author:   ted.schlicke  $
% $Date:   Feb 02 2018 11:40:10  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
    help Mike.readDfs0File
    return
end

options=struct;
options.info=false;
options.parameters=[];
options.locations=[];
options=checkArguments(options,varargin);

if ~exist(dfs0FileName,'file')
    error('File ''%s'' not found',dfs0FileName)
end

if ~stringFinder(dfs0FileName,'.dfs0','type','end','output','bool')
    error('File ''%s'' not a .dfs0 file',dfs0FileName)
end

timeScaleFactor=24*60*60;

if ~stringFinder(dfs0FileName,'.dfs0','type','end','output','bool')
    error('Not a dfs0 file! Aborting...')
end
%
NETaddAssembly('DHI.Generic.MikeZero.DFS.dll');
import DHI.Generic.MikeZero.DFS.*;
import DHI.Generic.MikeZero.DFS.dfs0.*;

dfs0File  = DfsFileFactory.DfsGenericOpen(dfs0FileName);

%% Dfs0 Items
NItems=dfs0File.ItemInfo.Count;
dfs0Items=arrayfun(@(i)char(dfs0File.ItemInfo.Item(i-1).Name),1:NItems,'Unif',0)';
% Data with location info seem to have format 'Location: Parameter'
% Data without location info have format 'Parameter'
% Split string into these two things:
% Locations & Parameters
str=regexp(dfs0Items,':','split');
str=vertcat(str{:});
if size(str,2)==1 % No location info
    dfs0Locations=[];
    dfs0Parameters=strtrim(unique(str));
else
    dfs0Locations=unique(str(:,1));
    dfs0Parameters=strtrim(unique(str(:,2)));
end
% Make 1st char lower-case
dfs0Parameters=cellfun(@(x)[lower(x(1)),x(2:end)],dfs0Parameters,'unif',0);
dfs0Items=cellfun(@(x)[lower(x(1)),x(2:end)],dfs0Items,'unif',0);

%% Read times and data for all items
% Use the Dfs0Util for bulk-reading all data and timesteps
dd = double(Dfs0Util.ReadDfs0DataDouble(dfs0File));

% Time stuff -
timeColumn = dd(:,1); % first column of data - seem to be in seconds past start time
% In dfs0 file generated from Mike -> File -> New, timestep seems to be
% minutes, even though dfs0File.FileInfo.TimeAxis.TimeUnit = 'eumUsec' ?!
t0=dfs0File.FileInfo.TimeAxis.StartDateTime;
y=t0.Year;
m=t0.Month;
d=t0.Day;
H=t0.Hour;
M=t0.Minute;
S=t0.Second;
vec=double([y,m,d,H,M,S]);
t0=datenum(vec);
% Use above parameters to generate timings of data (in datenums)
dateTime=t0+timeColumn/timeScaleFactor;

%% METADATA: Prepare struct with basic dfs0 information
metaData=struct('fileName',dfs0FileName,'fileTitle',char(dfs0File.FileInfo.FileTitle));
metaData.locations=dfs0Locations;
metaData.parameters=dfs0Parameters;
metaData.items=dfs0Items;

dfs0Output=struct('metaData',metaData);

% Does user just want info about file?
if options.info
    fprintf('Returning info only\n')
    return
end

%% Check data requirements
if ~isempty(options.parameters)
    fprintf('Checking parameters...\n')
    kParam=stringFinder(dfs0Items,options.parameters,'output','index','type','or','ignorecase',1);
    if isempty(kParam)
        fprintf('Valid parameters:\n')
        disp(dfs0Parameters)
        error('Specified parameter not found; please select one of the above')
    end
else
    kParam=1:NItems; % All items
end

if ~isempty(options.locations)
    fprintf('Checking locations...\n')
    kLoc=stringFinder(dfs0Items,options.locations,'output','index','type','or','ignorecase',1);
    if isempty(kLoc)
        if ~isempty(dfs0Locations)
            fprintf('Valid locations:\n')
            disp(dfs0Locations)
        else
            fprintf('No location information available; please don''t filter by location\n')
        end
        error('Specified location not found')
    end
else
    kLoc=1:NItems;
end

itemIndices2Extract=intersect(kParam,kLoc);

pointData=dd(:,2:end);

%% Load each location and parameter's data into struct array
NItems=length(itemIndices2Extract);
dataStructArray=cell(NItems,1); % Convert to struct array after loop
for itemIndex=1:NItems
    ithItemIndex=itemIndices2Extract(itemIndex);
    %    fprintf('Loading item %d of %d...\n',i,NItems)
    item=dfs0File.ItemInfo.Item(ithItemIndex-1);
    name=char(item.Name);  % name given by MIKE, for example "Shuna_Castle_Bay_2000: Surface elevation"
    %        unit=char(item.Quantity.Unit);
    splitName=strsplit(name,':'); % First bit of name (Observation point name)
    if length(splitName)==1
        locationName=[];
        paramName=strtrim(splitName{1});
    else
        locationName=strtrim(splitName{1});
        paramName=strtrim(splitName{end});
    end
    paramName(1)=lower(paramName(1));
    quantity=item.Quantity;
    unitDescription=char(quantity.UnitDescription);
    unitAbbreviation=char(quantity.UnitAbbreviation);
    %    fieldName=genvarname(sprintf('%s(%s)',paramName,unitAbbreviation));
    s=struct('dfs0Name',name,'location',locationName,'parameter',paramName,'dateTime',dateTime,'value',pointData(:,ithItemIndex),'unitDescription',unitDescription,'unitAbbreviation',unitAbbreviation);
    dataStructArray{ithItemIndex}=s;
end
dataStructArray=vertcat(dataStructArray{:}); %cell -> struct array
dfs0Output.data=dataStructArray;

dfs0File.Close()

