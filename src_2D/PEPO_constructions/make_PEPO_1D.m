function obj = make_PEPO_1D(obj)
    d = obj.dim;
    ln_prefact = obj.nf;
     obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
    obj.bounds = [1];
    obj.boundary_vect(obj.bounds) = 1;

    for n = 2:obj.order

        [map, ~] = create_map(1:n, obj.numopts);

        if mod(n, 2) == 1
            m = (n - 1) / 2;
            pattern = {[m, 0, m, 0]};
        else
            obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^(n)];
            obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^(n)];
            obj.current_max_index = n / 2;

            m = n / 2;
            pattern = {[m - 1, 0, m, 0], [m, 0, m - 1, 0]};
        end

        [obj, target, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

        e1 = calculate_error(obj, 1:n + 1, obj.numopts);
        %e2 = svds(target, 1);
        fprintf("n=%d residual = %.4e   d_nf %.4e  \n", n, e1,  exp(ln_prefact - obj.nf));

        if rank_x == 0%ineffective step, truncate

            if mod(n, 2) ~= 1
                obj.virtual_level_sizes_horiz = obj.virtual_level_sizes_horiz(1:end - 1);
                obj.virtual_level_sizes_vert = obj.virtual_level_sizes_vert(1:end - 1);
                obj.current_max_index = obj.current_max_index - 1;
            end

            obj.max_index = obj.current_max_index;
            break;
        end
    end
end
