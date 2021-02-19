function obj = scale_PEPO_pattern(obj, pattern, fact)
    n = numel(pattern);

    for ii = 1:n
        pat = pattern{ii} + 1;
        obj.PEPO_cell{pat(1), pat(2), pat(3), pat(4)} = obj.PEPO_cell{pat(1), pat(2), pat(3), pat(4)} .* (fact.^(1 / n));
    end
end
