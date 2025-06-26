function [dfsOutput ] = readDfsFile( dfsFileName,varargin )
% Wrapper function for readDfs0File / readDfsuFile
%
% INPUT: Name of dfs0/dfsu file(s)
% 
% Optional inputs:
% varargin - passed to readDfs0File/readDfsuFile (see their help sections
% for details)
%
% OUTPUT:
% struct (array) as returned by read functions above
% 

dfsFileName=cellstr(dfsFileName);
N=length(dfsFileName);
dfsOutput=cell(N,1);

for fileIndex=1:N
    f=dfsFileName{fileIndex};
    [~,~,ext]=fileparts(f);
    switch ext
        case '.dfs0'
            s=Mike.readDfs0File(f,varargin);
        case '.dfsu'
            s=Mike.readDfsuFile(f,varargin);
        otherwise
            error('Invalid valid')
    end
    dfsOutput{fileIndex}=s;
end

try
    dfsOutput=vertcat(dfsOutput{:});
catch
    % Fails if mixed dfs0/dfsu input
end