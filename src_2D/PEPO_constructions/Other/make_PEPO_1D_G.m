function [obj, err_code] = make_PEPO_1D_type_E(obj)
    err_code = 0;

    obj.complex = 1;

    d = obj.dim;
    ln_prefact = obj.nf;
    obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
    obj.bounds = [1];
    obj.boundary_vect(obj.bounds) = 1;

    m_i = floor(obj.order / 2);
    prime_level = m_i;

    obj.virtual_level_sizes_horiz = obj.dim.^(2 * [0:m_i, 1:m_i]);
    obj.virtual_level_sizes_vert = [1];

    assert(obj.max_index) = numel(obj.virtual_level_sizes_horiz);

    %obj.max_index = numel(obj.virtual_level_sizes_horiz);
    obj.current_max_index = obj.max_index;

    %obj.PEPO_cell{1, 1, 2 + prime_level, 1} = obj.PEPO_cell{1, 1, 2, 1};

    for n = 2:obj.order

        [map, ~] = create_map(1:n, obj.numopts);

        if mod(n, 2) == 1
            m = (n - 1) / 2;

            pattern = {[m + prime_level, 0, m + prime_level, 0]};
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

        else

            %obj.current_max_index = n / 2;

            m = n / 2;

            obj.PEPO_cell{m, 1, 1 + m + prime_level, 1} = 0.5 * reshape(eye(d^(2 * m)), [d, d, d^(2 * m - 2), 1, d^(2 * m), 1]) / exp(obj.nf);
            obj = assign_perm(obj, [m - 1, 0, m + prime_level, 0]);

            pattern = {[m - 1, 0, m, 0], [m, 0, m - 1, 0]};
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

        end

        %e1 = calculate_error(obj, 1:n, obj.numopts);
        %e2 = svds(target, 1);
        %fprintf("n=%d residual = %.4e  tar %.4e d_nf %.4e  \n", n, e1, e2, exp(ln_prefact - obj.nf));

    end
end
