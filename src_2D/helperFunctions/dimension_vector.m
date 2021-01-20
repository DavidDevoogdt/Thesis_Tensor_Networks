function p = dimension_vector(d, n, leftright)
    %helper function to create a 1xn vector  [ left,d,d,..,d,right]
    %if left/right are not supplied/0, this is omitted

    if nargin < 3
        p = zeros(1, n);
        p = p + d;
        return

    else

        p = zeros(1, n + 2);
        p = p + d;
        p(1) = leftright(1);
        p(end) = leftright(2);
    end

end
