function obj = make_PEPO_1D_double(obj)
    d = obj.dim;
    ln_prefact = obj.nf;

    
    
    m_i = floor(obj.order/2);
    prime_level = m_i ;

    obj.virtual_level_sizes_horiz = obj.dim.^(2 * [0:m_i, 1:m_i]);
    obj.virtual_level_sizes_vert = [1];
    
    assert(obj.max_index)= numel(obj.virtual_level_sizes_horiz);
    
    %obj.max_index = numel(obj.virtual_level_sizes_horiz);
    obj.current_max_index = obj.max_index;

    obj.PEPO_cell{1, 1, 2 + prime_level, 1} = obj.PEPO_cell{1, 1, 2, 1};

    for n = 2:obj.order

        [map, ~] = create_map(1:n, obj.numopts);

        if mod(n, 2) == 1
            m = (n - 1) / 2;
            if n == 3
                pattern = {[m, 0, m + prime_level, 0], [m + prime_level, 0, 0, 0]};
            else
                pattern = {[m, 0, m + prime_level, 0], [m + prime_level, 0, m + prime_level - 1, 0]};
            end
        else

            %obj.current_max_index = n / 2;

            m = n / 2;

            if n == 2
                pattern = {[0, 0, m, 0], [m, 0, 0, 0]};
            else
                pattern = {[m - 1, 0, m, 0], [m, 0, m - 1 + prime_level, 0]};
            end
        end

        [obj, target, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

        %e1 = calculate_error(obj, 1:n, obj.numopts);
        %e2 = svds(target, 1);
        %fprintf("n=%d residual = %.4e  tar %.4e d_nf %.4e  \n", n, e1, e2, exp(ln_prefact - obj.nf));

    end
end
