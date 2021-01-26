function [H_exp, lnprefact] = H_exp(obj, map, lnprefact, improve_ln_pre_factor)

    H_mat = H_matrix(obj, map);

    d = obj.dim;

    if nargin < 3
        lnprefact = obj.nf;
    end

    if nargin < 4
        improve_ln_pre_factor = false;
    end

    d_log_nf = 1;

    while abs(d_log_nf) > 1e-10
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
