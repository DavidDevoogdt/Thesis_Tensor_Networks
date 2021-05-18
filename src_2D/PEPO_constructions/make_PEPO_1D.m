function [obj, err_code] = make_PEPO_1D(obj)

    err_code = 0;

    d = obj.dim;
    ln_prefact = obj.nf;
    obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
    obj.bounds = [1];
    obj.boundary_vect(obj.bounds) = 1;

    for n = 2:obj.copts.order

        [map, ~] = create_map(1:n, obj.numopts);

        if mod(n, 2) == 1
            m = (n - 1) / 2;
            pattern = {[m, 0, m, 0]};
        else
            obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^(n)];
            obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^(n)];

            m = n / 2;
            pattern = {[m - 1, 0, m, 0], [m, 0, m - 1, 0]};
        end

        [obj, ln_prefact, err] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

        if obj.testing == 1
            e1 = calculate_error(obj, 1:n + 1, obj.numopts);
            fprintf("n=%d residual = %.4e   d_nf %.4e  \n", n, e1, exp(ln_prefact - obj.nf));
        end

    end
end
