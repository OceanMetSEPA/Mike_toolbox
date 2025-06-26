function op=massPerTimeStep(sourceInfo)
% Calculate model mass per timestep

NSources=height(sourceInfo);

op=nan(NSources,1);
p='\\asb-fp-mod01\AMMU\MarineModelling\Projects\Dervaig\Model\Run1_Mesh4.2.1_HDWind.m21fm - Result Files\';
for sourceIndex=1:NSources
    isource=sourceInfo(sourceIndex,:);
    %    fprintf('Processing source %d of %d\n',sourceIndex,NSources)
    %    disp(isource)
    % Key parameters for scaling:
    asFlux=isource.asFlux;
    perParticle=isource.perParticle;
    NParticles=isource.number_of_particles_per_timestep;
    particleMass=isource.particleMass;
    dt=isource.timeStep;

    dfs0File=char(isource.dfs0File);
    if ~isempty(dfs0File)
        dfs0File=char(GetFullPath(fullfile(p,dfs0File)));
        if ~isfile(dfs0File)
            error('File ''%s'' not found!',dfs0File)
        end
        dfs0=Mike.readDfs0File(dfs0File);
        data=dfs0.data;
%        t=data.dateTime;
        particlesPerTimeStep=data.value;
    else % Not dfs0 input
        % Generate appropriate time sequence
%        t=t0+(ind0:dind:ind1)'*dtHD/spd;
%        particlesPerTimeStep=repmat(NParticles,size(t));
        particlesPerTimeStep=NParticles;
    end
    % Now we have number of particles per timestep, either from dfs0 file
    % or based on m21fm file settings.
    % Scale these to appropriate mass:
    if ~perParticle % Don't worry about the number of particles, just if there is one:
        particlesPerTimeStep=particlesPerTimeStep>0;
    end
    massPerTimeStep=particlesPerTimeStep*particleMass;
    if asFlux
        massPerTimeStep=massPerTimeStep*dt;
    end
    op(sourceIndex)=massPerTimeStep;

end

% and we're done!