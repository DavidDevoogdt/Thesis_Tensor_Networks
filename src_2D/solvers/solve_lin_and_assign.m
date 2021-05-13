function [obj, target_site, res_target, ln_prefact_out, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact, opts)

    d = obj.dim;

    if ~isfield(opts, 'loop_dim')
        opts.loop_dim = -1;
    end
    if ~isfield(opts, 'loop')
        opts.loop = 0;
    end

    if ~( isfield(opts, 'target_site')  && isfield(opts, 'all_con_cells') && isfield(opts, 'pat_cells'))

        [all_con_cells, pat_cells] = get_valid_contractions(obj, map, struct('max_index', obj.current_max_index, 'pattern', {pattern}));
        
        [target, ln_prefact_out] = H_exp(obj, map, ln_prefact, true);
        target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(d^2, map.N));
        
        con_cells = all_con_cells(pat_cells);
        
        if map.N > 0 %do matrix contraction and substract non contracted cells

            obj = cell2matrix(obj);

            target = permute(reshape(target_site, dimension_vector(d, 2 * map.N)), [1:2:2 * map.N, 2:2:2 * map.N]);

            target = target - contract_network(obj, map, struct('matrix', 1, 'lnprefact', ln_prefact_out));
            target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(obj.dim^2, map.N));

            target_site = -contract_con_cells(obj, map, ln_prefact_out, -target_site, con_cells);

        else
            target_site = contract_con_cells(obj, map, ln_prefact_out, target_site, all_con_cells(~pat_cells));
        end
    else
        
        target_site = opts.target_site;
        ln_prefact_out = ln_prefact;
       
        
        all_con_cells = opts.all_con_cells;
        pat_cells = opts.pat_cells;
        
        con_cells = all_con_cells(pat_cells);
        
    end

    if isfield(opts, 'sub_I')
        target_site = target_site - 1 * reshape(eye(d^map.N), dimension_vector(d^2, map.N)) / exp(ln_prefact_out * map.N);

    end

    mul_factor = exp(ln_prefact_out - obj.nf);

    

    [x_cell, res_target, rank_x, res_con] = solve_lin(obj, pattern, map, con_cells, target_site, ln_prefact_out, opts.loop_dim, opts.loop);

    if rank_x ~= 0
        for i = 1:size(x_cell, 2)
            obj.PEPO_cell{pattern{i}(1) + 1, pattern{i}(2) + 1, pattern{i}(3) + 1, pattern{i}(4) + 1} = x_cell{i} * mul_factor;

            if obj.testing == 2
                fprintf("%.4e ", norm(reshape(x_cell{i} * mul_factor, [], 1), 2));
            end
        end
    end

    res_target = ipermute(reshape(res_target, dimension_vector(d, 2 * map.N)), site_ordering_permute(map.N));
    res_target = reshape(res_target, [d^map.N, d^map.N]);

end
