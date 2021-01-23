function [obj, target, res_target, ln_prefact_out, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact, loop_dim)
    if nargin <5
        loop_dim = -1;
    end

    d = obj.dim;

    [target, ln_prefact_out] = H_exp(obj, map, ln_prefact, true);
    target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(d^2, map.N));
   

    mul_factor = exp(ln_prefact_out - obj.nf);

    con_cells = get_valid_contractions(obj, map, struct('max_index', obj.current_max_index, 'pattern', {pattern}));

    [x_cell, res_target, rank_x,res_con] = solve_lin(obj, pattern, map, con_cells, target_site, ln_prefact_out,loop_dim);

    
    if rank_x ~= 0
        for i = 1:size(x_cell, 2)
            obj.PEPO_cell{pattern{i}(1) + 1, pattern{i}(2) + 1, pattern{i}(3) + 1, pattern{i}(4) + 1} = x_cell{i} * mul_factor;
        end
    end
    
    res_target = ipermute(reshape(res_target, dimension_vector(d, 2 * map.N)), site_ordering_permute(map.N));
     
    if obj.testing ==1
        temp_list_1 = fetch_PEPO_cells(obj, map, res_con{1}{1}, ln_prefact_out);
        A1 = ncon(temp_list_1, map.leg_list);
        diff = A1-res_target;
    end

    target = reshape(target, [d^map.N, d^map.N]);
    res_target = reshape(res_target, [d^map.N, d^map.N]);

    % if obj.testing == 1
    %     err = calculate_error( obj,1:n, obj.numopts);
    %     print(err);
    % end
end
