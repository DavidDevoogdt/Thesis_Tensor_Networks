function [obj, err_code] = make_PEPO_1D(obj)

    err_code = 0;

    d = obj.dim;
    ln_prefact = obj.nf;
    obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
    obj.bounds = [1];
    obj.boundary_vect(obj.bounds) = 1;
    
    m_i = floor(obj.copts.order / 2);
    prime_level = m_i;

    obj.virtual_level_sizes_horiz = obj.dim.^(2 * [0:m_i, 1:m_i]);
    

    for n = 2:obj.copts.order

        [map, ~] = create_map(1:n, obj.numopts);

        if mod(n, 2) == 1
            m = (n - 1) / 2;
            pattern = {[m, 0, m, 0]};
            [obj, ln_prefact, err] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct() );

            
        else
            
            %sz = min(d^n, obj.copts.max_bond_dim);
            
            %obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, sz];
            %obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, sz];

            m = n / 2;
            pattern = {[m - 1, 0, m, 0], [m, 0, m - 1, 0]};
            [obj, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct('remove_S',1) );

            pattern = {[m - 1, 0, m+prime_level, 0], [m+prime_level, 0, m - 1, 0]};
            [obj, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct );

            
            
        end

       
        if obj.testing == 1
            e1 = calculate_error(obj, 1:n + 1, obj.numopts);
            fprintf("n=%d residual = %.4e   d_nf %.4e  \n", n, e1, exp(ln_prefact - obj.nf));
        end

    end
end
