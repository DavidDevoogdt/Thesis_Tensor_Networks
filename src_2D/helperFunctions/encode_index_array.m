function y = encode_index_array(n, len, max_num)
    % this takes a single number and splits it into the composite indices. The
    i = 1;
    y = zeros(len, 1);

    while n ~= 0
        [n, r] = readOne(n, max_num + 1);
        y(i) = r;
        i = i + 1;
    end
end

function [s, r] = readOne(s, d)
    r = rem(s, d);
    s = (s - r) / d;
end
