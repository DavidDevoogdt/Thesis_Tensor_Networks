function p = site_ordering_permute(n, alt)

    if nargin < 2
        alt = 0;
    end

    % changes from |left i1 i2 ... j1 j2.. right> to |left i1 j1 i2 j2 ...
    % right>
    p = zeros(2 * n, 1);
    %p(1)=1;
    %p(2*n+2)=2*n+2;
    if alt == 0

        for i = 1:n
            p(2 * i - 1) = i;
        end

        for i = 1:n
            p(2 * i) = n + i;
        end
    else
        for i = 1:n
            p(2 * i) = i;
        end

        for i = 1:n
            p(2 * i - 1) = n + i;
        end
    end
end
