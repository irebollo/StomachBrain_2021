
function [z,mu,sigma] = nanzscore(x,flag,dim)

% [] is a special case for std and mean, just handle it out here.
if isequal(x,[]), z = x; return; end

if nargin < 2
    flag = 0;
end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

z = NaN(size(x));
inan = ~isnan(x);
x = x(inan);
[ztmp,mu,sigma] = zscore(x,flag,dim);
z(inan) = ztmp;