function [err] = calculate_error(obj, nummap, opts, matrix, non_trace_num)

    if isfield(nummap, 'map')
        map = nummap;
        b_map = map;
    else
        [map, b_map] = create_map(nummap, opts);
    end

    d = obj.dim;

    cn_opts.max_index = obj.max_index;

    if nargin < 4
        cn_opts.matrix = 0;
    else
        cn_opts.matrix = matrix;
    end

    if cn_opts.matrix == 1
        obj = cell2matrix(obj);
    end

    if nargin == 5
        cn_opts.trace = true;
        cn_opts.non_trace_num = non_trace_num;
    else
        cn_opts.trace = false;
    end

    Contraction = contract_network(obj, map, cn_opts);

    H_matrix = H_exp(obj, map, obj.nf);

    if cn_opts.trace == true
        b = permute(H_matrix, site_ordering_permute(map.N2));
        convect = [reshape([1:non_trace_num - 1; 1:non_trace_num - 1], [], 1)', -1, -2, reshape([non_trace_num:map.N2 - 1; non_trace_num:map.N2 - 1], [], 1)'];
        b = ncon({b}, {convect});

        a = reshape(permute(Contraction, site_ordering_permute(map.N2)), [d, d]);

        err = norm((b - a), 2) / norm(b, 2);
        %prefact =
    else
        err = calculate_error_core(H_matrix-Contraction , H_matrix);
        
%         b = reshape(H_matrix, [d^(map.N2), d^(map.N2)]);
%         a = reshape(Contraction, [d^(map.N2), d^(map.N2)]);
% 
%         p = 2;
% 
%         [~, S1, ~] = svds(a - b, 30);
% 
%         sum_1 = (sum(diag(S1).^p))^(1 / p);
% 
%         [~, S2, ~] = svds(b, 30);
% 
%         sum_2 = (sum(diag(S2).^p))^(1 / p);
% 
%         %prefact = exp(obj.nf * map.N) * (sum_2);
% 
%         % for real values, multiply both with obj.nf^(map.N)*trace_a =
%         % prefact^N
% 
%         err = sum_1 / sum_2;
    end
end
