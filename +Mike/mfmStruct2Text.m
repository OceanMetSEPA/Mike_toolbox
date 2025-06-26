function lines = mfmStruct2Text(ms, numberOfSpaces)
% Recursive function to return cell array of strings representing struct contents

if nargin < 2
    numberOfSpaces = 0;
    clc
end

fn = fieldnames(ms);
lines = {};  % Always initialize as a cell array

% Detect top-level call to process only uppercase fields
initialCall = (numberOfSpaces == 0);
if initialCall
    lines=cellfun(@(x)[x,sprintf('\r')],ms.header,'unif',0);
    engineName = ms.EngineName;
    lines{end+1,1} = sprintf('[%s]\r', engineName);

    % Only process upper case fields
    fn = fn(strcmp(fn, upper(fn)));
end

Nf = length(fn);
NSpaces = numberOfSpaces;

for fieldIndex = 1:Nf
    fni = fn{fieldIndex};
    spaceString = repmat(' ', 1, numberOfSpaces + 3);
    val = ms.(fni);

    if isstruct(val)
        lines{end+1,1} = sprintf('%s[%s]\r', spaceString, fni);

        % Recursive call
        subLines = Mike.mfmStruct2Text(val, NSpaces + 3);

        % FORCE subLines to be a cell array of strings
        if ischar(subLines)
            subLines = {subLines};
        elseif isnumeric(subLines)
            subLines = {num2str(subLines)};
        elseif ~iscell(subLines)
            subLines = {tdisp(subLines)};
        end
        try
            lines = [lines; subLines(:)];  % Ensure column
        catch err
            disp(err)
            fprintf('Lines: %s\n',tdisp(size(lines)))
            disp(lines)
            underline
            fprintf('subLines: %s\n',tdisp(size(lines)))
            disp(subLines)
            error('error!')
        end
        lines{end+1,1} = sprintf('%sEndSect  // %s\r', spaceString, fni);
        lines{end+1,1} = sprintf('\r');
    else
        if isnumeric(val)
            valString = strjoin(arrayfun(@(x) sprintf('%.15g', x), val(:)', 'UniformOutput', false), ', ');
        else
            valString = tdisp(val);
        end
        lines{end+1,1} = sprintf('%s%s = %s\r', spaceString, fni, valString);
    end
end

if initialCall
    lines{end+1,1} = sprintf('EndSect  // %s\r', engineName);
end

% Force output to always be a cell array of char
if ~iscell(lines)
    lines = {lines};
end
end
