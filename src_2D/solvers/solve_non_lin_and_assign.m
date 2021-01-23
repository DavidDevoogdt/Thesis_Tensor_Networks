function [obj, ln_prefact] = solve_non_lin_and_assign(obj, map, pattern, ln_prefact,loop_dim,start_size)

    if nargin<6
        start_size = 0.1;
    end

    d = obj.dim;

    [target, ln_prefact] = H_exp(obj, map, ln_prefact, true);
    target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(d^2, map.N));

    mul_factor = exp(ln_prefact - obj.nf);

    %guess initial value
    for i = 1:size(pattern, 2)
        pattern_s = [d,d,  (pattern{i} ~=0)*loop_dim +  (pattern{i} ==0)];
        obj.PEPO_cell{pattern{i}(1) + 1, pattern{i}(2) + 1, pattern{i}(3) + 1, pattern{i}(4) + 1} = rand(pattern_s)*start_size;
    end


    con_cells = get_valid_contractions(obj, map, struct('max_index', obj.current_max_index, 'pattern', {pattern}));

    x_cell = solve_non_lin(obj, pattern, {map}, {target_site}, {con_cells}, struct(), ln_prefact  )


    for i = 1:size(pattern, 2)
        obj.PEPO_cell{pattern{i}(1) + 1, pattern{i}(2) + 1, pattern{i}(3) + 1, pattern{i}(4) + 1} = x_cell{i} * mul_factor;
    end

end
