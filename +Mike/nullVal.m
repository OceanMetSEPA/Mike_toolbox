function nullValue=nullVal()
% This works if called from command line, but not from script:
% nameOfThisScript = matlab.desktop.editor.getActiveFilename;
% p=fileparts(nameOfThisScript);
% nullFile=fullfile(p,'mikeNullValue.mat');
% nullValue=importdata(nullFile);
% stk = dbstack;
% filepath = which(stk(1).file);
% 
% for i=1:length(stk)
%     disp(stk(i))
%     underline
% end

tmp=which('Mike.nullVal');
tmp=strrep(tmp,'nullVal.m','mikeNullValue.mat');
nullValue=importdata(tmp);

