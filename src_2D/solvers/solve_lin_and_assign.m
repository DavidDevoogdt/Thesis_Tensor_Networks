function [obj, ln_prefact_out, err] = solve_lin_and_assign(obj, map, pattern, ln_prefact, opts)
    %linear solver one provide single map, patterns should be one or 2
    %neighbouring patters. Otherwise use solve_sequential_lin_and_assing or
    %solve_non_lin_and_assign. opts can be struct
    %Full inverse is legacy

    d = obj.dim;

    split_opts = struct;
    
    if isfield(opts, 'svd_split_dim')
        split_opts.svd_split_dim = opts.svd_split_dim;
    end
    if isfield(opts, 'remove_S')
        split_opts.remove_S = opts.remove_S;
    end
    
    inv_opts = struct;
    
    if ~isfield(opts, 'full_inverse')
        inv_opts.full_inverse = 0;
    else
        inv_opts.full_inverse = opts.full_inverse ;
    end

   
    
    
    %contract con_cells for map, if not provided as argument
    if ~(isfield(opts, 'target_site') && isfield(opts, 'all_con_cells') && isfield(opts, 'pat_cells'))

        [all_con_cells, pat_cells] = get_valid_contractions(obj, map, struct('pattern', {pattern}));
        [target, ln_prefact_out] = H_exp(obj, map, ln_prefact, true);
        con_cells = all_con_cells(pat_cells);

        if nargout == 3
            target_orig = target;
        end

        %do matrix contraction and substract non contracted cells if problem is large enough. Saves time for large N
        if map.N > 5

            obj = cell2matrix(obj);

            target = target - contract_network(obj, map, struct('matrix', 1, 'lnprefact', ln_prefact_out));
            target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(obj.dim^2, map.N));

            target_site = -contract_con_cells(obj, map, ln_prefact_out, -target_site, con_cells);

        else
            target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(d^2, map.N));
            target_site = contract_con_cells(obj, map, ln_prefact_out, target_site, all_con_cells(~pat_cells));
        end
    else
        target_site = opts.target_site;
        if nargout == 3
            target_orig = permute(reshape(target_site, dimension_vector(d, 2 * map.N)), [1:2:2 * map.N, 2:2:2 * map.N]);
        end
        ln_prefact_out = ln_prefact;
        all_con_cells = opts.all_con_cells;
        pat_cells = opts.pat_cells;
        con_cells = all_con_cells(pat_cells);
    end

    mul_factor = exp(ln_prefact_out - obj.nf);
    
    inv_opts.inv_eps = obj.copts.inv_eps*mul_factor;
    
    %if isfield(opts, 'mul_factor')
        %split_opts.mul_factor_n = mul_factor^map.N;
        %split_opts.sigma =  obj.copts.err_tol;
    %end

    %actually compute inverse
    [x_cell, res_target, res_con] = solve_lin(obj, pattern, map, con_cells, target_site, ln_prefact_out,split_opts , inv_opts);

    %assign cells
    for i = 1:size(x_cell, 2)
        obj.PEPO_cell{pattern{i}(1) + 1, pattern{i}(2) + 1, pattern{i}(3) + 1, pattern{i}(4) + 1} = x_cell{i} * mul_factor;

        if obj.testing == 2
            fprintf("%.4e ", norm(reshape(x_cell{i} * mul_factor, [], 1), 2));
        end
    end

    %compute error efficiently
    if nargout == 3
        [~, res_target] = optimize_con_cells(obj, {map}, {res_con}, {}, {res_target}, ln_prefact_out);
        res_target = permute(reshape(res_target{1}, dimension_vector(d, 2 * map.N)), [1:2:2 * map.N, 2:2:2 * map.N]);
        err = calculate_error_core(res_target, target_orig);
    end

end
