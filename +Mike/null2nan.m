function op = null2nan(ip)
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
if isstruct(ip)
    ip=replaceValInStruct(ip,oldVal,newVal);
elseif isnumeric(ip)
    ip(ip==oldVal)=newVal;
else
    error('Input should be struct or numeric')
end
op=ip;
