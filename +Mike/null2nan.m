function s = null2nan(s)
%NULL2NAN Recursively replace numeric occurrences of nullVal with NaN in a struct
%
% INPUT:
%   s       - input struct (possibly nested, with struct arrays or cell arrays)
%
% OUTPUT:
%   s       - struct with nullVal replaced by NaN in numeric fields
%
% Example:
%   nullVal=Mike.nullVal;
%   s = struct('a',nullVal, 'b', [1,nullVal], 's', struct('x',nullVal));
%   sOut = null2nan(s);

oldVal=Mike.nullVal;
newVal=nan;

s=replaceValInStruct(s,oldVal,newVal);