function [b] = skipnans(a,idx)

b = NaN(size(idx));
for i = 1:numel(idx)
    if ~isnan(idx(i))
        b(i) = a(idx(i));
    end
end