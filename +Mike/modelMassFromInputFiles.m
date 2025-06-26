function op=modelMassFromInputFiles(mfmStruct,plotit)
% Calculate model mass (kg) from struct containing info about input file (.m21fm, dfs0...)
%
% INPUT: mfmStruct - struct as returned by Mike.mfm2struct()
%
% OUTPUT: struct with fields:
% Name - model run name
% dateTime - datenums of timesteps
% mass - timeseries of cumulative mass
% massPerTimeStep - mass added each timestep
% timeFilter - used in conv function to calculate mass
% 


if length(mfmStruct)>1
    op=arrayfun(@(i)Mike.modelMassFromInputFiles(mfmStruct(i)),1:length(mfmStruct),'unif',0)';
    try
        op=vertcat(op{:});
    catch
    end
    return
end

if ~isstruct(mfmStruct)
    error('Need struct input')
end
% Maybe settings struct we require is field of another struct? 
if isfield(mfmStruct,'m21fmStruct')
    mfmStruct=mfmStruct.m21fmStruct;
end

if nargin<2
    plotit=false;
end

% Time info:
spd=24*60*60;
dtHD=mfmStruct.TIME.time_step_interval;
ind0=mfmStruct.HYDRODYNAMIC_MODULE.DECOUPLING.first_time_step;
ind1=mfmStruct.HYDRODYNAMIC_MODULE.DECOUPLING.last_time_step;
%fprintf('Time steps %d to %d, interval = %d\n',ind0,ind1,dtHD)
dind=1;
t0=datenum(mfmStruct.TIME.start_time);
% Mike fm input filename and path:
f=mfmStruct.FileName;
p=fileparts(f);

% Find source info:
sourceInfo=Mike.sourceInfo(mfmStruct);
NSources=height(sourceInfo);

if plotit
    cm=rainbow(NSources);
    prepareFigure('close',1)
    h=nan(NSources,1);
end

op=cell(NSources,1);

for sourceIndex=1:NSources
    isource=sourceInfo(sourceIndex,:);
    %    fprintf('Processing source %d of %d\n',sourceIndex,NSources)
    %    disp(isource)
    % Key parameters for scaling:
    asFlux=isource.asFlux;
    perParticle=isource.perParticle;
    NParticles=isource.number_of_particles_per_timestep;
    particleMass=isource.particleMass;

    dfs0File=char(isource.dfs0File);
    if ~isempty(dfs0File)
        dfs0File=char(GetFullPath(fullfile(p,dfs0File)));
        if ~isfile(dfs0File)
            error('File ''%s'' not found!',dfs0File)
        end
        dfs0=Mike.readDfs0File(dfs0File);
        data=dfs0.data;
        t=data.dateTime;
        particlesPerTimeStep=data.value;
    else % Not dfs0 input
        % Generate appropriate time sequence
        t=t0+(ind0:dind:ind1)'*dtHD/spd;
        particlesPerTimeStep=repmat(NParticles,size(t));
    end
    % Now we have number of particles per timestep, either from dfs0 file
    % or based on m21fm file settings.
    % Scale these to appropriate mass:
    if ~perParticle
        particlesPerTimeStep=particlesPerTimeStep>0;
    end
    massPerTimeStep=particlesPerTimeStep*particleMass;
    if asFlux
        massPerTimeStep=massPerTimeStep*dtHD;
    end
    % CHECK
    % check=repmat(Mike.massPerTimeStep(isource),length(t),1);
    % if ~isequal(check,massPerTimeStep)
    %     check
    %     massPerTimeStep
    %     error('MISMATCH!')
    % end

    % Now we have mass of particles per timestep.
    % So we can calculate mass within model, factoring in decay.
    % Prepare time filter:
    decayFactorPerSecond=sourceInfo.decay(sourceIndex); % /s
    decayFactorPerDay=decayFactorPerSecond*spd; % /days
    timeFilter=exp(-decayFactorPerDay*(t-t(1)));

    % % Calculate mass in model
    % Convolution works for non-uniform input:
    mass=conv(massPerTimeStep,timeFilter);
    assignin('base','massPerTimeStep',massPerTimeStep)
    assignin('base','timeFilter',timeFilter)
    assignin('base','massFromConv',mass)
    % % Crop to appropriate length (Has length 2Nt-1 for equal length
    % inputs)
    % Also, conv applies input at start of timestep. Want 1st one to be
    % zero, so offset:
    mass=[0;mass(1:length(t)-1)];
    assignin('base','mass',mass)

    % Theoretical mass with exponential decay (CONSTANT MASS PER TIMESTEP):
    massPerSecond=mean(massPerTimeStep)/dtHD;
    k=decayFactorPerSecond;
    tsec=(t-t(1))*spd;
    mexp=massPerSecond/k*(1-exp(-k*tsec));

    % Store mass vs time for this source
    op{sourceIndex}=struct('Name',char(isource.Name),'dateTime',t,'mass',mass,'massPerTimeStep',massPerTimeStep,'timeFilter',timeFilter);
    if plotit
        h(sourceIndex)=plot(t,mass,'color',cm(sourceIndex,:),'displayname',char(isource.Name));
        plot(t,mexp,'--b','linewidth',2)
    end
end

if plotit
    legHandle=legend(h);
    set(legHandle,'interpreter','none')
    datetimeAxis
    zoom on
end

% Bundle up output:
op=vertcat(op{:});

% and we're done!