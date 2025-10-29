function op=nan2null(ip)
% Replace nan values in struct with Mike's null value (1.00000001800251e-35)
%
% INPUT:
% s - struct / array
%
% OUTPUT:
% s - struct / array with nan values in numeric fields replaced with null value
%
% Note - Mike uses its null value for e.g. particle positions outwith
% grid. This is useful when calling Mike.meshIndex for example (which
% returns nan in such cases), since matlab's tsearchn function (key
% component of meshIndex function) doesn't support nan values. But
% sometimes a nan value is preferable e.g. if we want to plot tracks -
% points outwith mesh aren't displaced.
%
% EXAMPLE:
% s=struct('val',[1:5,nan,6:10])
% Mike.nan2null(s) % returns struct with fields  val: [1 2 3 4 5 1.00000001800251e-35 6 7 8 9 10]

if nargin==0
    help Mike.nan2null
    return
end

newVal=Mike.nullVal;

if isstruct(ip)
    ip=replaceValInStruct(ip,nan,newVal);
elseif isnumeric(ip)
    ip(isnan(ip))=newVal;
else
    error('Input should be struct or numeric')
end
op=ip;
