function op=modelMassFromOutputFiles(dfsu,meshStruct)
% % Calculate model mass (kg) from output
%
% INPUTS:
% dfsu - struct as returned by Mike.readDfsuFile
% meshStruct - used by sussed2Total function if required
%
% OUTPUT:
% Name - model run name
% dfsuField - field of dfsu struct used for calculating mass
% dateTime - datenums of model output
% mass - timeseries of total mass in model domain
%
% NB total calculated incorrectly in Mike. Use sussed2Total function
% to regenerate but requires total water depth field, not always output...
%

dfsuParametersToProcess=stringFinder(fieldnames(dfsu),'sussed');
if isempty(dfsuParametersToProcess)
    dfsu=Mike.sussed2Total(dfsu,meshStruct);
end
NParameters2Process=length(dfsuParametersToProcess);

op=cell(NParameters2Process,1);

t=dfsu.dateTime;
for sourceIndex=1:NParameters2Process
    fni=dfsuParametersToProcess{sourceIndex};
    % Get source name from parameter:
    sourceName=fni;
%    sourceName=strrep(sourceName,'total','');
    sourceName=strrep(sourceName,'Mass','');
    sourceName=strrep(sourceName,'0x2C','');
    sourceName=strrep(sourceName,'sussedsum','');
    fprintf('%d and fni = %s; source name = %s\n',sourceIndex,fni,sourceName)
    vals=dfsu.(fni);
    % Convert total concentration to mass
    vals=sum(vals);
    vals=vals/1e9; % Convert to kg
    vals=vals(:);
    op{sourceIndex}=struct('Name',sourceName,'dfsuField',fni,'dateTime',t,'mass',vals);

end

op=vertcat(op{:});