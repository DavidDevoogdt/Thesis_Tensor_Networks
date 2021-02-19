function [err, prefact] = calculate_error(obj, nummap, opts,matrix)
    if nargin<4
        matrix = 0;
    end

    [map, b_map] = create_map(nummap, opts);

    d = obj.dim;

    H_matrix = H_exp(obj, map, obj.nf);

    if matrix==1
        Contraction = contract_network(obj, map, struct('max_index', obj.max_index, "matrix", 1));
    else
        Contraction = contract_network(obj, b_map, struct('max_index', obj.max_index));
    end

    b = reshape(H_matrix, [d^(map.N2), d^(map.N2)]);
    a = reshape(Contraction, [d^(map.N2), d^(map.N2)]);

    p = 2;

    [~, S1, ~] = svds(a - b, 30);

    sum_1 = (sum(diag(S1).^p))^(1 / p);

    [~, S2, ~] = svds(b, 30);

    sum_2 = (sum(diag(S2).^p))^(1 / p);

    prefact = exp(obj.nf * map.N) * (sum_2);

    % for real values, multiply both with obj.nf^(map.N)*trace_a =
    % prefact^N

    err = sum_1 / sum_2;
end
