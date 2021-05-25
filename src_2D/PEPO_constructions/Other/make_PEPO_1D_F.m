function [obj, err_code] = make_PEPO_1D(obj)

    err_code = 0;

    d = obj.dim;
    ln_prefact = obj.nf;
    obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
    obj.bounds = [1];
    obj.boundary_vect(obj.bounds) = 1;

    obj.complex = true;

    rot_180 = {{[3, 4, 1, 2]}};
    nl_opts = struct('Gradient', true, 'Display', 'none', 'maxit', 30);

    %[obj, ln_prefact] = solve_non_lin_and_assign(obj, { create_map( [1,1] , struct( 'h_cyclic',1  )) }, { [0,0,0,0]  }, ln_prefact, struct('Gradient', false, 'Display', 'none', 'maxit', 20));

    %obj.PEPO_cell{1, 1, 1, 1} = reshape((expm(obj.H_1_tensor)) / exp(obj.nf), [d, d, 1, 1, 1, 1]);

    function [obj, ln_prefact, err] = add_lin(obj, pattern, ln_prefact, nl)
        if nargin < 4
            nl = 0;
        end

        [map1, ~] = create_map(make_cross(pattern));

        if nl == 0
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, {pattern}, ln_prefact);
        else
            [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map1}, {pattern}, ln_prefact, nl_opts, rot_180);
        end
        obj = assign_perm(obj, pattern);

        err = calculate_error(obj, map1, [], 1);

        if obj.testing == 1
            if err > obj.err_tol
                disp(err);
            end
        end
    end

    for n = 2:obj.order

        if mod(n, 2) == 1
            m = (n - 1) / 2;

            [obj, ln_prefact, ~] = add_lin(obj, [m, m, 0, 0], ln_prefact);

        else

            obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^(n)];
            obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^(n)];
            obj.current_max_index = n / 2;

            m = n / 2;

            [obj, ln_prefact, ~] = add_lin(obj, [m - 1, 0, m, 0], ln_prefact, 1);
        end

        err02 = calculate_error(obj, 1:10, obj.cycleopts, 1);

    end
end
