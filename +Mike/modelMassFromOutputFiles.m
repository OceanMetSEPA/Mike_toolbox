function op=modelMassFromOutputFiles(dfsu,meshStruct,varargin)
% % Calculate model mass (kg) from MIKE output file
%
% INPUTS:
% dfsu - struct as returned by Mike.readDfsuFile
% meshStruct - used by sussed2Total function if required
%
% Optional Inputs:
%   massField ['sussed'] - calculate mass of field(s) containing this
%   scaleFactor [1e9] - scale mass in µg
%
% OUTPUT:
% struct with fields:
% *) Name - model run name
% *) dfsuField - field of dfsu struct used for calculating mass
% *) dateTime - datenums of model output
% *) mass - timeseries of total mass in model domain
%
% NB total calculated incorrectly in Mike. Use sussed2Total function
% to regenerate but requires total water depth field, not always output...
%
% Alternatively, from particle tracking we may calculate concs ourselves
% (e.g. sealice stuff). In which case we need to calculate mass by
% multiplying concs by area of cells

options=struct;
options.scaleFactor=1e9;
options.massField='sussed';
options=checkArguments(options,varargin);

fprintf('Finding mass within dfsu file...\n')

dfsuParametersToProcess=stringFinder(fieldnames(dfsu),options.massField);
if isempty(dfsuParametersToProcess)
    fprintf('Combining suspended/sedimented concentrations...\n')
    dfsu=Mike.sussed2Total(dfsu,meshStruct);
end
dfsuParametersToProcess=stringFinder(fieldnames(dfsu),options.massField);
NParameters2Process=length(dfsuParametersToProcess);

fprintf('Calculating mass of these parameters:\n')
disp(dfsuParametersToProcess)

op=cell(NParameters2Process,1);

t=dfsu.dateTime;
for sourceIndex=1:NParameters2Process
    fni=dfsuParametersToProcess{sourceIndex};
    % Get source name from parameter:
    sourceName=fni;
    sourceName=strrep(sourceName,'Mass','');
    sourceName=strrep(sourceName,'0x2C','');
    sourceName=strrep(sourceName,'sussedsum','');
    fprintf('%d and fni = %s; source name = %s\n',sourceIndex,fni,sourceName)
    vals=dfsu.(fni);
    % Convert total concentration to mass
    vals=sum(vals);
    vals=vals/options.scaleFactor; % Convert to kg
    vals=vals(:);
    op{sourceIndex}=struct('Name',sourceName,'dfsuField',fni,'dateTime',t,'mass',vals);
end

op=vertcat(op{:});
