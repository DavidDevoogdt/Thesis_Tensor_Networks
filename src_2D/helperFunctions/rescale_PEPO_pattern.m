function obj = rescale_PEPO_pattern(obj, pattern, fact, num)
    %assumes every pattern only occurs once

    if nargin < 3
        fact = 1.5;
        num = 0;
    end

    n = numel(pattern);

    sizes = zeros(n, 1);

    for i = 1:n
        pat = pattern{i} + 1;
        sizes(i) = max(reshape(abs(obj.PEPO_cell{pat(1), pat(2), pat(3), pat(4)}), [], 1));
    end

    [min_argvalue, argmin] = min(sizes);
    [max_argvalue, argmax] = max(sizes);

    sqrt = (min_argvalue * max_argvalue)^0.5;

    %make max smaller
    max_pat = pattern{argmax} + 1;
    min_pat = pattern{argmin} + 1;

    obj.PEPO_cell{max_pat(1), max_pat(2), max_pat(3), max_pat(4)} = obj.PEPO_cell{max_pat(1), max_pat(2), max_pat(3), max_pat(4)} .* (sqrt / max_argvalue);
    obj.PEPO_cell{min_pat(1), min_pat(2), min_pat(3), min_pat(4)} = obj.PEPO_cell{min_pat(1), min_pat(2), min_pat(3), min_pat(4)} .* (sqrt / min_argvalue);

    sizes(argmin) = sqrt;
    sizes(argmax) = sqrt;

    if max(sizes) / min(sizes) > fact
        obj = rescale_PEPO_pattern(obj, pattern, fact, num + 1);
    end

    if num == 0
        for i = 1:n
            pat = pattern{i} + 1;
            sizes(i) = max(reshape(abs(obj.PEPO_cell{pat(1), pat(2), pat(3), pat(4)}), [], 1));
        end
        if obj.testing == 1
            sizes
        end
    end
end
