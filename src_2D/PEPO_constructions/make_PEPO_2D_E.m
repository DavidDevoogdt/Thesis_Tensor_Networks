function [obj, error_code] = make_PEPO_2D_A(obj)

    load_save = 0;

    function [obj, ln_prefact, err] = add_lin(obj, pattern, ln_prefact, ldim)
        if nargin < 4
            ldim = -1;
        end

        if numel(pattern) == 2
            %only assign this pattern
            [map1, ~] = create_map(make_cross(pattern{1}));
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, pattern, ln_prefact, struct("loop_dim", ldim));

        else
            perms = get_perm(pattern{1});

            for i = 1:size(perms, 1)
                dperms = perms{i, 1};
                for j = 1:size(dperms, 1)
                    ppp = dperms(j, :);
                    if obj.testing == 1
                        disp(ppp)
                    end
                    if isempty(obj.PEPO_cell{ppp(1) + 1, ppp(2) + 1, ppp(3) + 1, ppp(4) + 1})
                        [map1, ~] = create_map(make_cross(ppp));
                        [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, {ppp}, ln_prefact, struct);
                    end
                end
            end
        end

        err = check_err(obj, map1);
        
        
    end

    function [obj, ln_prefact] = add_loop(obj, arr, ln_prefact)

        perms = get_perm([arr, b, a], [0, 0, 1, 1]);

        for i = 1:size(perms, 1)
            dperms = perms{i, 1};
            for j = 1:size(dperms, 1)
                ppp = dperms(j, :)
                if isempty(obj.PEPO_cell{ppp(1) + 1, ppp(2) + 1, ppp(3) + 1, ppp(4) + 1})

                    [map1, ~] = create_map(make_cross_loop(ppp, ppp > L));

                    %[obj, ln_prefact] = solve_non_lin_and_assign(obj, {map1}, {ppp}, ln_prefact, nl_opts);

                    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, {ppp}, ln_prefact, struct);
                end
            end
        end

        check_err(obj, map1)

        %obj = cell2matrix(obj); max(reshape( obj.PEPO_matrix ,[],1))

    end

    if load_save == 0

        %setup and size definitions
        d = obj.dim;
        ln_prefact = obj.nf;

        error_code = 0;

        rot_180 = {{[3, 4, 1, 2]}};
        rot_90 = {[2, 3, 4, 1], ...
                [3, 4, 1, 2], ...
                [4, 1, 2, 3]};

        %obj.testing = 1

        if obj.testing == 1
            nl_opts = struct('Gradient', true, 'Display', 'iter', 'maxit', 40);
        else
            nl_opts = struct('Gradient', true, 'Display', 'none', 'maxit', 40);
        end

        obj.current_max_index = 1;
        obj.virtual_level_sizes_horiz = [1];
        obj.virtual_level_sizes_vert = [1];

        obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
        obj.bounds = [1];
        obj.boundary_vect(obj.bounds) = 1;

        obj.complex = true;

        L = obj.L;
        max_bond_dim = obj.max_bond_dim;

        for n = 1:L
            sz = min(d^(2 * n), max_bond_dim);

            obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, sz];
            obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, sz];
        end

        %% function definitions for adding linear and loop blocks

        %%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% Blocks 1D chain %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if obj.testing == 1
            fprintf('linear blocks')
        end

        %blocks (same as 1D chain), checks whether cyclic error improves with extra added blocks
        obj.PEPO_cell{1, 1, 1, 1} = reshape((expm(obj.H_1_tensor)) / exp(obj.nf), [d, d, 1, 1, 1, 1]);

        %[obj, ln_prefact] = solve_non_lin_and_assign(obj, { create_map( [1, 1] , struct) }, { [0,0,0,0]  }, ln_prefact, struct('Gradient', false, 'Display', 'iter', 'maxit', 20))

        for n = 1:L

            if obj.testing == 1
                fprintf('linear: %d', n)
            end
            %%% n-1--|--n--|--n1- block

            if n ~= 1
                err01 = calculate_error(obj, 1:(2 * n + 2), obj.cycleopts, 1);
            else
                err01 = 1;
            end

            [obj, ln_prefact, ~] = add_lin(obj, {[n - 1, 0, n, 0], [n, 0, n - 1, 0]}, ln_prefact, obj.virtual_level_sizes_horiz(n + 1));

            %rotate
            obj.PEPO_cell{1, n, 1, n + 1} = permute(obj.PEPO_cell{n, 1, n + 1, 1}, [1, 2, 6, 3, 4, 5]);
            obj.PEPO_cell{1, n + 1, 1, n} = permute(obj.PEPO_cell{n + 1, 1, n, 1}, [1, 2, 6, 3, 4, 5]);

            %do all other variants
            if n ~= 1
                [obj, ln_prefact, ~] = add_lin(obj, {[n - 1, 0, n, 0]}, ln_prefact);
            end

            err02 = calculate_error(obj, 1:(2 * n + 2), obj.cycleopts, 1);

            if err02 > err01
                warning('not converging')
                error_code = 1;

                obj.PEPO_cell{n + 1, n + 1, 1, 1} = [];
                obj = assign_perm(obj, [n, n, 0, 0]);
            end

            %%% n--|--n block

            if n < 3

                [obj, ln_prefact, ~] = add_lin(obj, {[n, n, 0, 0]}, ln_prefact);
                err03 = calculate_error(obj, 1:(2 * n + 2), obj.cycleopts, 1);

                if err03 > err02
                    warning('not converging')
                    error_code = 1;
                    obj.PEPO_cell{n + 1, n + 1, 1, 1} = [];
                    obj = assign_perm(obj, [n, n, 0, 0]);
                end
            end
        end

        %%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% blocks with 3/4 legs %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if obj.testing == 1
            fprintf('non linear blocks')
        end

        for n = 1:L
            switch n
                case 1
                    [obj, ln_prefact, ~] = add_lin(obj, {[1, 1, 1, 0]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[1, 1, 1, 1]}, ln_prefact);

                case 2
                    [obj, ln_prefact, ~] = add_lin(obj, {[2, 1, 1, 0]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[2, 1, 1, 1]}, ln_prefact);

                    [obj, ln_prefact, ~] = add_lin(obj, {[2, 2, 1, 0]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[2, 2, 1, 1]}, ln_prefact);

                    [obj, ln_prefact, ~] = add_lin(obj, {[2, 2, 2, 0]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[2, 2, 2, 1]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[2, 2, 2, 2]}, ln_prefact);
                case 3
                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 1, 1, 0]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 1, 1, 1]}, ln_prefact);

                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 2, 1, 0]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 2, 1, 1]}, ln_prefact);

                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 2, 2, 0]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 2, 2, 1]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 2, 2, 2]}, ln_prefact);

                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 3, 1, 0]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 3, 1, 1]}, ln_prefact);

                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 3, 2, 0]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 3, 2, 1]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 3, 2, 2]}, ln_prefact);

                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 3, 3, 0]}, ln_prefact);
                    [obj, ln_prefact, ~] = add_lin(obj, {[3, 3, 3, 1]}, ln_prefact);

                    %[obj, ln_prefact, ~] = add_lin(obj, {[3, 3, 3, 2]}, ln_prefact); %4*60s
                    %[obj, ln_prefact, ~] = add_lin(obj, {[3, 3, 3, 3]}, ln_prefact); %starts swapping, n
                otherwise
                    error('not impl')
            end
        end
    else
        load('temp_save.mat')
    end
    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% LOOPS %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    a = numel(obj.virtual_level_sizes_horiz);
    b = a + 1;
    c = b + 1;
    e = c + 1;

    obj.current_max_index = e;
    obj.max_index = e;

    dd = [6, 8, 8];

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, dd];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, dd];

    single_loop = 1;
    single_extensions = 0;
    offset_loop = 0;
    double_extension = 0;
    double_loop = 0;
    L_loop = 0;

    if single_loop == 1
        if obj.testing == 1
            fprintf('simple loop\n')
        end

        [map, ~] = create_map([
                            1, 1;
                            1, 1; ], struct);

        alpha_pattern_perm = {rot_90};

        pattern = {[a, a, 0, 0]};

        [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);
        check_err(obj, map);
    end

    if single_extensions == 1

        if obj.testing == 1
            fprintf('single_ext\n')
        end

        for k = 0:3

            [m, pattern] = rotate([
                                0, 0, 0, 0;
                                1, 1, 1, 0;
                                0, 1, 1, 0; ], {[0, b, a, 0], [c, 0, 0, a], [1, 0, c, b]}, k);

            [map, ~] = create_map(m, struct);
            obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            obj = rescale_PEPO_pattern(obj, pattern);

            %[obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts);

            check_err(obj, map);

            [m, pattern] = rotate([
                                0, 1, 0, 0;
                                0, 1, 1, 0;
                                0, 1, 1, 0; ], {[0, 1, c, b]}, k);

            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

            [m, pattern] = rotate([
                                0, 0, 0, 0, 0;
                                1, 1, 1, 1, 0;
                                0, 0, 1, 1, 0; ], {[2, 0, c, b]}, k);

            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

            [m, pattern] = rotate([
                                0, 1, 0, 0;
                                0, 1, 0, 0;
                                0, 1, 1, 0;
                                0, 1, 1, 0; ], {[0, 2, c, b]}, k);

            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

            [m, pattern] = rotate([
                                0, 1, 0, 0;
                                1, 1, 1, 0;
                                0, 1, 1, 0; ], {[1, 1, c, b]}, k);

            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

            [m, pattern] = rotate([
                                0, 1, 0, 0;
                                0, 1, 0, 0;
                                1, 1, 1, 0;
                                0, 1, 1, 0; ], {[1, 2, c, b]}, k);

            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

            [m, pattern] = rotate([
                                0, 0, 0, 0, 0;
                                0, 0, 1, 0, 0;
                                1, 1, 1, 1, 0;
                                0, 0, 1, 1, 0; ], {[2, 1, c, b]}, k);

            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

            [m, pattern] = rotate([
                                0, 0, 1, 0, 0;
                                0, 0, 1, 0, 0;
                                1, 1, 1, 1, 0;
                                0, 0, 1, 1, 0; ], {[2, 2, c, b]}, k);

            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

            %loops with 1 or more legs
            %
            %     for n = 1:L
            %         switch n
            %             case 1
            %                 [obj, ln_prefact] = add_loop(obj, [1, 0], ln_prefact);
            %                 %[obj, ln_prefact] = add_loop(obj, [1, 1], ln_prefact);
            %             case 2
            %                 %[obj, ln_prefact] = add_loop(obj, [2, 0], ln_prefact);
            %                 %[obj, ln_prefact] = add_loop(obj, [2, 1], ln_prefact);
            %                 %[obj, ln_prefact] = add_loop(obj, [2, 2], ln_prefact);
            %             case 3
            %                 [obj, ln_prefact] = add_loop(obj, [3, 0], ln_prefact);
            %                 [obj, ln_prefact] = add_loop(obj, [3, 1], ln_prefact);
            %                 [obj, ln_prefact] = add_loop(obj, [3, 2], ln_prefact);
            %                 [obj, ln_prefact] = add_loop(obj, [3, 3], ln_prefact);
            %             otherwise
            %                 error('not impl')
            %         end
            %     end

        end

    end

    if offset_loop == 1

        if obj.testing == 1
            fprintf('offset loops\n')
        end
        for k = 0:1
            [m, pattern] = rotate([
                                1, 1, 0, 0;
                                1, 1, 1, 0;
                                0, 1, 1, 0; ], {[c, b, c, b]}, k);

            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);
            %[obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts);

            check_err(obj, map);

        end
    end

    if double_extension == 1
        if obj.testing == 1
            fprintf('double extension\n')
        end
        for k = 1:4
            [m, pattern] = rotate([0, 0, 0, 0;
                                1, 1, 1, 1;
                                0, 1, 1, 0; ], {[1, 0, b, a], [b, 0, 1, a]}, k);

            [map, ~] = create_map(m, struct);

            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct('loop', 1, 'loop_dim', obj.virtual_level_sizes_horiz(b + 1)));
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern(1), ln_prefact, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern(2), ln_prefact, struct);

            check_err(obj, map);

            [m, pattern] = rotate([0, 0, 1, 0;
                                1, 1, 1, 0;
                                0, 1, 1, 0; ], {[b, 1, 0, a]}, k);
            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

            [m, pattern] = rotate([0, 0, 1, 0;
                                1, 1, 1, 1;
                                0, 1, 1, 0; ], {[b, 1, 1, a]}, k);
            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

            [m, pattern] = rotate([0, 1, 0, 0;
                                0, 1, 1, 1;
                                0, 1, 1, 0; ], {[0, 1, b, a]}, k);
            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

            [m, pattern] = rotate([0, 1, 0, 0;
                                1, 1, 1, 1;
                                0, 1, 1, 0; ], {[1, 1, b, a]}, k);
            [map, ~] = create_map(m, struct);
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

            check_err(obj, map);

        end
    end

    if double_loop == 1

        if obj.testing == 1
            fprintf('double loop\n')
        end

        for k=0:1
             [m, pattern] = rotate([
                                1, 1, 1 ;
                                1, 1, 1; ], {
                                %[0, 0, b, a], [b, 0, c, a],[c,0,0,a],[0,a,b,0],[b,a,c,0],[c,a,0,0]
                                [a,0,a,b],[a,b,a,0]
                                }, k);
            
            [map, ~] = create_map(m, struct);

            obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, struct('display', obj.testing));
            obj = rescale_PEPO_pattern(obj, pattern); %

            check_err(obj, map);
        end
      
    end

    if L_loop == 1

        if obj.testing == 1
            fprintf('L shape \n')
        end

        for k = 1:4
            [m, pattern] = rotate([1, 1, 1;
                                1, 1, 1;
                                1, 1, 0; ], {[b, b, a, a]}, k);

            [map, ~] = create_map(m, struct);
            tic
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);
            toc
            check_err(obj, map);
        end
    end

    %     calculate_error(obj, [
    %                     0, 0, 0, 0;
    %                     0, 1, 1, 0;
    %                     0, 1, 1, 0;
    %                     0, 0, 0, 0; ], [], 1);

end

function err = check_err(obj, map1)
    err = calculate_error(obj, map1, [], 1);

    if obj.testing == 1
        if err > obj.err_tol
            disp(err);
        end
    end

    obj = cell2matrix(obj);
    msize = max(  abs( reshape(obj.PEPO_matrix,[],1)   )  );
    if msize>2
        fprintf('size %.4e',msize);
    end
    
end


function [m, p] = rotate(m, p, k)
    m = rot90(m, k);
    for l = 1:numel(p)
        p{l} = circshift(p{l}, -k);
    end
end
