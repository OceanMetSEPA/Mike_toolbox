function op=sourceInfo(varargin)
% Extract key info from mfmStruct
% (which is generated from m21fm file by Mike.mfm2struct(fileName);
%
% INPUT:
% mfmStruct - struct returned my Mike.mfm2struct, or .m21fm (etc) input file
%
% OUTPUT:
% table with fields:
%   *) Name - source name
%   *) coordinates - location of source (x,y,z)
%   *) ClassIndex - number of corresponding class
%   *) number_of_particles_per_timestep
%   *) particleMass 
%   *) unit - of mass ('g' or 'kg' for now)
%   *) perParticle - is mass divided between particles?
%   *) asFlux - is mass released per second? (otherwise per timestep)
%   *) dfs0NParticles - file containing time-varying release
%   *) dfs0Mass - file containing time-varying mass
%   *) decay (per second) - rate at which particles decay
%   *) timeStep (s) - model timestep (used for calculating mass elsewhere)
%

mfmStruct=[varargin{:}];

if ischar(mfmStruct)
    if ~isfile(mfmStruct)
        error('file ''%s'' not fount',mfmStruct)
    end
    mfmStruct=Mike.mfm2struct(mfmStruct);
end
N=length(mfmStruct);
if N>1
    op=arrayfun(@(i)Mike.sourceInfo(mfmStruct(i)),1:N,'unif',0);
    op=vertcat(op{:});
    return
end

NSources=mfmStruct.PARTICLE_TRACKING_MODULE.SOURCES.number_of_sources;
if ischar(NSources) % mfm2struct didn't convert to numeric
    mfmStruct=numstruct(mfmStruct); % so do it now
end
s=mfmStruct.PARTICLE_TRACKING_MODULE;
NSources=s.SOURCES.number_of_sources;


op=cell(NSources,1);
for i=1:NSources
    fni=sprintf('SOURCE_%d',i);
    source=s.SOURCES.(fni);
    if ~source.include
        continue
    end
    sourceName=strrep(source.Name,'''','');
    xyz=source.coordinates;
    % check classes for this source. DHI stores info about all classes for
    % each source. The one we're after in this instance seems to be the one
    % where constant_value >0
    % 01/08/2024 - this fails for ArdowBurn run - 2 classes have
    % constant_value = 1000???
    % k=find(arrayfun(@(i)s.SOURCES.(fni).(sprintf('CLASS_%d',i)).constant_value,1:NClasses,'unif',1));
    k=i;
    switch length(k)
        case 0
            % No class for this source- maybe it's not included
            continue
        case 1
        otherwise
            error('too many classes for this source')
    end
    % If we're here, we've got one class
    classStr=sprintf('CLASS_%d',k);
    iclass=s.SOURCES.(fni).(classStr);
    % Prepare struct with key info for this source
    istruct=struct;
    istruct.Name=sourceName;
    istruct.coordinates=xyz;
    istruct.ClassIndex=k;
    istruct.number_of_particles_per_timestep=iclass.number_of_particles_per_timestep;
    istruct.particleMass=iclass.constant_value;
    istruct.unit=Mike.eumCode2Unit(mfmStruct.PARTICLE_TRACKING_MODULE.CLASSES.(classStr).EUM_unit);
    istruct.perParticle=~iclass.type_particle;
    istruct.asFlux=logical(iclass.type_value);
    if iclass.format_particle
        istruct.dfs0NParticles=strrep(iclass.file_name_particle,'|','');
        istruct.dfs0Mass=strrep(iclass.file_name,'|','');
    else
        istruct.dfs0NParticles='';
        istruct.dfs0Mass='';
    end
    % Decay
    if mfmStruct.PARTICLE_TRACKING_MODULE.DECAY.(classStr).type
        istruct.decay=mfmStruct.PARTICLE_TRACKING_MODULE.DECAY.(classStr).constant_value;
    else
        istruct.decay=0;
    end
    % Timestep (useful to calculating mass)
    istruct.timeStep=mfmStruct.TIME.time_step_interval;
    op{i}=istruct;
end

op=vertcat(op{:});

try % struct2table doesn't like tables with one row
    op=struct2table(op);
catch
    op=struct2table(op,'AsArray',1);
end
% switch nargout
%     case 0
%         disp(sourceInfo)
%     case 1
%         varargout{1}=sourceInfo;
%     otherwise
%         error('too many outputs')
% end
end
