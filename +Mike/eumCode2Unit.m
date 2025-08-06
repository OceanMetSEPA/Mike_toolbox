function unit=eumCode2Unit(x)
% Find EUM description from code
% 
% INPUT:
% x - scalar integer
%
% OUTPUT:
% unit - name of unit
%
% Mike input definition files (e.g. .m21fm) use codes rather than
% descriptions for the units of mass in particle tracking classes. Useful
% to know what these are!
%
% Might want to expand this to get more info 
% e.g. unit type (particles are of type 'Component Mass', code 100039).  
% For now, just get masses so we can scale concentrations properly
%
if ~isnumeric(x)
    error('Numeric input required')
end
if ~isscalar(x)
    error('scalar input required')
end

% Define unit names, codes, and scale factors relative to gram
% Codes found from MikeIO in python with help from ChatGPT
% (might want co
%unitNames = {'kilogram'; 'gram'; 'milligram'; 'microgram'};
unitNames = {'kg'; 'g'; 'mg'; 'µg'};
unitCodes = [1200; 1201; 1202; 1203];
% Create table with the additional scale factor column
%unitsTable = table(unitNames, unitCodes,'VariableNames', {'Unit', 'Code'});

k=x==unitCodes;
if ~any(k)
    error('No matching codes found!')
end

unit=unitNames{k};