function [H_exp, lnprefact] = H_exp(obj, map, lnprefact, improve_ln_pre_factor)
    %calculate the tensor exponential throug matrix diagonalisation for the given map.
    %improve_ln_pre_factor finds an optimal normalistion factor, at the cost of doing the diagonalistion multiple times

    H_mat = H_matrix(obj, map);

    d = obj.dim;

    if nargin < 3
        lnprefact = obj.nf;
    end

    if nargin < 4
        improve_ln_pre_factor = false;
    end

    d_log_nf = 1;

    while abs(d_log_nf) > 0.1
        H_matrix_2 = H_mat - eye(d^(map.N)) * map.N * lnprefact;

        H_expo = expm(H_matrix_2);
        d_log_nf = log(svds(H_expo, 1)) / map.N;

        lnprefact = lnprefact + d_log_nf;

        if ~improve_ln_pre_factor
            break
        end
    end

    lnprefact = lnprefact - d_log_nf;

    H_exp = reshape(H_expo, dimension_vector(d, 2 * map.N));
end
