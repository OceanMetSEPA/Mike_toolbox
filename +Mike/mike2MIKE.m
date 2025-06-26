function s=mike2MIKE(s)
% Convert structs created from dfs0/dfsu files created by Mike into format
% created by MIKE, for backward compatability
%
% Specifically:
% 1) Convert field/parameter names from camelCase (capitalise first letter)
% 2) Bundle struct array containing data (one struct per parameter) into
% single struct
%
% INPUT:
% s - struct created by Mike.readDfs(0/u)File
%
% s - struct as created by MIKE.readDfs(0/u)File
%
% EXAMPLE
% f='\\asb-fp-mod01\AMMU\MarineModelling\Projects\EastSkye\Model\test3d_Alt7_p2.m3fm - Result Files\SCQN1_17042020_WL.dfs0';
% old=MIKE.readDfs0File(f);
% new=Mike.readDfs0File(f);
% isequal(old,new); % 0
% redo=Mike.mike2MIKE(new);
% isequal(old,redo) % 1

% Fix fieldnames:
s=unCamelCase(s);

% Sort data:
if isfield(s,'MetaData')
    if ~isempty(s.MetaData.Locations)
        if isfield(s,'Data')
            s.Data=dataArray2Struct(s.Data);
        end
    end
end
% Regenerate MetaData.Items - fixes minor inconsistency in camelCase naming
if isfield(s,'MetaData')
    % Final fiddle to make sure items consistent
    locs=s.MetaData.Locations;
    if ~isempty(locs)
        N=length(s.MetaData.Items);
        for i=1:N
            item=s.MetaData.Items{i};
            bits=strsplit(item,':');
            bit1=bits{1};
            loc=locs{strcmpi(bit1,locs)};
            bits{1}=loc;
            item=sprintf('%s:%s',bits{1},bits{2});
            s.MetaData.Items{i}=item;
        end
    end
end
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sub-functions - used to be in separate files but probably don't need to
% call them directly..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s=unCamelCase(s,varargin)
% Convert struct fieldnames and char values from camelCase (i.e. capitalise
% first letter)
%
% For units, capitalise fieldname but not value
%
% Also, remove ',' from fieldnames if present

options=struct;
%options.ignore={'UnitDescription','UnitAbbreviation','Items','Location','Locations'};
options.ignore={'UnitDescription','UnitAbbreviation','Location','Locations'};
options=checkArguments(options,varargin);

fn=fieldnames(s);
nfn=cellfun(@(x)[upper(x(1)),x(2:end)],fn,'unif',0);

% Recursive bit for fields that themselves are structs
NStructs=length(s);
ca=cell(NStructs,1);
for structIndex=1:NStructs
    si=s(structIndex);
    NFields=length(nfn);
    for fieldIndex=1:NFields
        fi=fn{fieldIndex};
        nfi=nfn{fieldIndex};
        nfi=strrep(nfi,'0x2C',''); % Remove any commas
        if strcmp(fi,nfi)
            continue
        end
        %        fprintf('Adding field %s and removing %s\n',nfi,fi)
        si.(nfi)=si.(fi);
        si=rmfield(si,fi);
        if ismember(nfi,options.ignore) % Don't convert units
            continue
        end
        val=si.(nfi);
        %        fprintf('Fiddling field %s\n',nfi)
        switch class(val)
            case 'char'
                val=capitalise(val);
            case 'cell'
                for cellIndex=1:length(val)
                    ival=val{cellIndex};
                    ival={capitalise(ival)};
                    val{cellIndex}=ival;
                end
                val=vertcat(val{:});
            case 'struct'
                val=unCamelCase(val);
        end
        si.(nfi)=val;
    end
    ca{structIndex}=si;
end
s=vertcat(ca{:});
end

function str=capitalise(str)
exceptions={'x','y','z'};
if ismember(str,exceptions)
    return
end
str=[upper(str(1)),str(2:end)];
end

function dataStruct=dataArray2Struct(dataArray)
% Convert data struct array created by Mike to format as returned by MIKE
%
% Mike returns struct array with Nl*Np structs, where:
% Nl = number of locations
% Np = number of parameters
%
% Want structs for each location to be bundled together

if length(dataArray)==1 % Only need fiddle below if struct array
    return
end

dataArray=unCamelCase(dataArray);

charLocs=true;
switch class([dataArray.Location])
    case 'double'
        locs=[dataArray.Location];
        charLocs=false;
    case 'char'
        locs={dataArray.Location};
    otherwise
        error('unsupported Location type %s',class([dataArray.Location]))
end

ulocs=unique(locs);
Nl=length(ulocs);
dataStruct=cell(Nl,1);
for locIndex=1:Nl
    if charLocs
        loc=ulocs{locIndex};
        k=strcmp(loc,locs);
    else
        loc=ulocs(locIndex);
        k=loc==locs;
    end
    locData=dataArray(k);
    s=struct('Location',locData(1).Location,'DateTime',locData(1).DateTime);
    Np=length(locData);
    for parameterIndex=1:Np
        idata=locData(parameterIndex);
        % Prepare fieldname based on parameter
        param=idata.Parameter;
        for letterIndex=2:length(param)
            if ~isletter(param(letterIndex-1))
                param(letterIndex)=upper(param(letterIndex));
            end
        end
        param=strrep(param,' ',''); % remove spaces
        param=sprintf('%s(%s)',param,idata.UnitAbbreviation); % add unit
        param=genvarname(param); % make suitable for fieldname
        value=idata.Value;
        s.(param)=value;
    end
    dataStruct{locIndex}=s;
end
dataStruct=vertcat(dataStruct{:});
end

