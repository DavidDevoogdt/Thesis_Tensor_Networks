function [U2, S2, V2] = expand_svd(U, S, V, n)
    dim = size(S, 1);
    U2 = zeros(dim, dim + n);
    V2 = zeros(dim, dim + n);
    S2 = zeros(dim + n);

    if dim < dim + n
        U2(1:dim, 1:dim) = U;
        V2(1:dim, 1:dim) = V;
        S2(1:dim, 1:dim) = S;
    else
        U2 = U(:, 1:dim + n);
        V2 = V(:, 1:dim + n);
        S2 = S(1:dim + n, 1:dim + n);

    end

end

