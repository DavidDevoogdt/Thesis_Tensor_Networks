function [obj, error_code] = make_PEPO_2D_sym(obj)

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
        nl_opts = struct('Gradient', true, 'Display', 'iter', 'maxit', 50);
    else
        nl_opts = struct('Gradient', true, 'Display', 'none', 'maxit', 100);
    end

    %obj.current_max_index = 1;
    obj.virtual_level_sizes_horiz = [1];
    obj.virtual_level_sizes_vert = [1];

    obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
    obj.bounds = [1];
    obj.boundary_vect(obj.bounds) = 1;

    %obj.complex = false;

    L = obj.copts.L;
    max_bond_dim = obj.copts.max_bond_dim;

    for n = 1:L
        sz = min(d^(2 * n), max_bond_dim);

        obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, sz];
        obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, sz];
    end

    %% function definitions for adding linear and loop blocks

    function [obj, ln_prefact, err] = add_lin(obj, pattern, ln_prefact, nl)
        if obj.testing == 1
            disp(pattern)
        end

        if nargin < 4
            nl = 0;
        end

        [map1, ~] = create_map(make_cross(pattern));

        if nl == 0
            [obj, ln_prefact, err] = solve_lin_and_assign(obj, map1, {pattern}, ln_prefact, struct);
            %[obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, {pattern}, ln_prefact, struct);
        elseif nl == 1
            %[obj, ln_prefact, err] = solve_sequential_lin_and_assign(obj, map1, {pattern}, ln_prefact, struct('display', 1, 'maxit', 150), @(x) assign_perm(x, pattern), 1);
            [obj, ln_prefact, err] = solve_sequential_lin_and_assign(obj, map1, {pattern}, ln_prefact, struct('display', obj.testing, 'maxit', 150), {rot_90});
        else
            [obj, ln_prefact, err] = solve_non_lin_and_assign(obj, {map1}, {pattern}, ln_prefact, nl_opts, rot_180);
        end
        obj = assign_perm(obj, pattern);

        %err = calculate_error(obj, map1, [], 1);

        obj.copts.inv_eps = max( obj.copts.inv_eps,  10* err);
        
        if obj.testing == 1
            if err > obj.copts.err_tol
                disp(err);
            end
        end
    end

    function [obj, ln_prefact] = add_loop(obj, arr, ln_prefact)
        pat = [arr, c, b];

        if obj.testing == 1
            disp(pat)
        end

        map1 = make_cross_loop(pat, [0, 0, 1, 1]);
        map1 = create_map(map1, struct);

        [obj, ln_prefact, err] = solve_lin_and_assign(obj, map1, {pat}, ln_prefact, struct);

        obj = assign_perm(obj, pat, [0, 0, 1, 1]);

        %err = calculate_error(obj, map1, [], 1);
        
        obj.copts.inv_eps = max( obj.copts.inv_eps,  10* err);

        if obj.testing == 1
            if err > obj.copts.err_tol
                disp(err);
            end
        end

        %obj = cell2matrix(obj); max(reshape( obj.PEPO_matrix ,[],1))

    end

    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Blocks 1D chain %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if obj.testing == 1
        fprintf('linear blocks')
    end
    %blocks (same as 1D chain), checks whether cyclic error improves with extra added blocks
    obj.PEPO_cell{1, 1, 1, 1} = reshape((expm(obj.H_1_tensor)) / exp(obj.nf), [d, d, 1, 1, 1, 1]);

    for n = 2:obj.copts.order

        err01 = calculate_error(obj, 1:n + 1, obj.cycleopts, 1);

        if mod(n, 2) == 0
            m = n / 2;

            if obj.testing == 1
                fprintf('linear: %d\n', m)
            end

            if n == 2
                [obj, ln_prefact, err02] = add_lin(obj, [m - 1, 0, m, 0], ln_prefact, 2);
            else
                [obj, ln_prefact, err02] = add_lin(obj, [m - 1, 0, m, 0], ln_prefact, 1);
            end

        else
            m = (n - 1) / 2;

            %%% n--|--n block
            [obj, ln_prefact, err02] = add_lin(obj, [m, m, 0, 0], ln_prefact);
        end

        %err02 = calculate_error(obj, 1:n + 1, obj.cycleopts, 1);

        if err02 > err01 && err02 > obj.copts.err_tol
            warning('not converging old %.4e new %.4e n %d', err01, err02, n)
            error_code = 1;
            return;

        end

        err01 = err02;
    end

    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% blocks with 3/4 legs %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if obj.testing == 1
        fprintf('non linear blocks \n')
    end

    for n = 3:obj.copts.order
        switch n
            case 3
                [obj, ln_prefact, ~] = add_lin(obj, [1, 1, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [1, 1, 1, 1], ln_prefact);

            case 4
                [obj, ln_prefact, ~] = add_lin(obj, [2, 1, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [2, 1, 1, 1], ln_prefact);
            case 5
                [obj, ln_prefact, ~] = add_lin(obj, [2, 2, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [2, 2, 1, 1], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [2, 2, 2, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [2, 2, 2, 1], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [2, 2, 2, 2], ln_prefact);
            case 6
                [obj, ln_prefact, ~] = add_lin(obj, [3, 1, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 1, 1, 1], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [3, 2, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 2, 1, 1], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [3, 2, 2, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 2, 2, 1], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 2, 2, 2], ln_prefact);

            case 7

                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 1, 1], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 2, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 2, 1], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 2, 2], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 3, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 3, 1], ln_prefact);
                

                %works but slow and high ram consumptio
                %[obj, ln_prefact, ~] = add_lin(obj, [3, 3, 3, 2], ln_prefact);
                %[obj,ln_prefact,~] = add_lin(obj, [3,3,3,3] ,ln_prefact);
                %takes to much memory
            otherwise
                error('not impl')
        end
    end

    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% LOOPS %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    a = numel(obj.virtual_level_sizes_horiz);
    b = a + 1;
    c = b + 1;
    e = c + 1;

    %obj.current_max_index = e;
    %obj.max_index = e;

    %dd = [8, 8, 6];
    dd = [6, 8 8, 10];

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, dd];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, dd];

    if obj.copts.do_loops == 1
        if obj.testing == 1
            fprintf('simple loop\n')
        end
        %simple loop
        [map, ~] = create_map([1, 2; 3, 4], obj.numopts);

        alpha_pattern_perm = {rot_90, rot_90, rot_90};
        pattern = {[a, a, 0, 0]};

        [obj, ln_prefact, error] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);
    end

    if obj.copts.loop_extension == 1

        if obj.testing == 1
            fprintf('loop extension')
        end

        pattern = {[c, 0, 0, a], ...
                    [0, b, a, 0]};

        obj = fill_rand(obj, pattern);

        [map, ~] = create_map([
                            1, 1, 1;
                            0, 1, 1; ], struct);
        pattern2 = {[1, 0, c, b]};
        [obj, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern2, ln_prefact, struct);

        obj = rescale_PEPO_pattern(obj, [pattern, pattern2]);

        obj = assign_perm(obj, [1, 0, c, b], [0, 0, 1, 1]);
        obj = assign_perm(obj, [c, 0, 0, a], [1, 0, 0, 1]);
        obj = assign_perm(obj, [0, b, a, 0], [0, 1, 1, 0]);

        err1 = calculate_error(obj, map, [], 1);
        
        obj.copts.inv_eps = max( obj.copts.inv_eps,  10* err1);

        if obj.testing == 1
            if err1 > obj.copts.err_tol
                disp(err1);
            end
        end

        %loops with 1 or more legs
        for n = 3:obj.copts.order
            switch n
                case 3
                    [obj, ln_prefact] = add_loop(obj, [1, 0], ln_prefact);
                    [obj, ln_prefact] = add_loop(obj, [1, 1], ln_prefact);
                case 4
                    [obj, ln_prefact] = add_loop(obj, [2, 0], ln_prefact);
                    [obj, ln_prefact] = add_loop(obj, [2, 1], ln_prefact);
                case 5
                    [obj, ln_prefact] = add_loop(obj, [2, 2], ln_prefact);
                case 6
                    [obj, ln_prefact] = add_loop(obj, [3, 0], ln_prefact);
                    [obj, ln_prefact] = add_loop(obj, [3, 1], ln_prefact);
                    [obj, ln_prefact] = add_loop(obj, [3, 2], ln_prefact);
                case 7
                    [obj, ln_prefact] = add_loop(obj, [3, 3], ln_prefact);
                otherwise
                    error('not impl')
            end
        end

    end

    if obj.copts.double_extension == 1

        [map, ~] = create_map([
                            0, 0, 1, 0;
                            1, 1, 1, 0;
                            0, 1, 1, 0; ], struct);
        %pattern = {[1, 0 ,e , b],[e, 1, 0, c]}
        %obj = solve_sequential_lin_and_assign(obj, map, pattern, ln_prefact, struct('display', 1,'maxit',5))
        %[obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts);

        pattern = {[1, 0, e, b], [e, 1, 0, c]};
        %[obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct( 'svd_split_dim', ) );
        obj = assign_perm(obj, pattern{1}, [1, 0, 0, 1]);

        err1 = calculate_error(obj, map, [], 1);

        if obj.testing == 1
            if err1 > obj.copts.err_tol
                disp(err1);
            end
        end

        %             [map, ~] = create_map([
        %                                 0, 0, 1, 0;
        %                                 1, 1, 1, 1;
        %                                 0, 1, 1, 0; ], struct);
        %             pattern = {[c, 1, 1, a]};
        %             [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);
        %             obj = assign_perm(obj, pattern{1}, [1, 0, 0, 1]);
        %
        %             err1 = calculate_error(obj, map, [], 1);
        %
        %             if obj.testing == 1
        %                 if err1 > obj.copts.err_tol
        %                     disp(err1);
        %                 end
        %             end

    end

    if obj.copts.offset_loops == 1

        if obj.testing == 1
            fprintf('offset loops\n')
        end

        [map, ~] = create_map([
                            1, 1, 0;
                            1, 1, 1;
                            0, 1, 1; ], struct);
        %pattern = {[b, a, c, b]};
        %pattern = {[b, a, a, c]};
        %pattern = {[c, b, c, b]};
        %pattern = {[c, b, a, a]};
        pattern = {[c, a, a, b]};

        [obj, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, struct);

        obj = scale_PEPO_pattern(obj, pattern, 0.5);
        obj = assign_perm(obj, pattern{1}, [1, 1, 1, 1]);

        err1 = calculate_error(obj, map, [], 1);

        if obj.testing == 1
            if err1 > obj.copts.err_tol
                disp(err1);
            end
        end
    end

    if obj.copts.double_loop == 1

        if obj.testing == 1
            fprintf('offset loops\n')
        end

        [map, ~] = create_map([
                            1, 1, 1;
                            1, 1, 1; ], struct);
        pattern = {[b, 0, a, c]};

        obj = solve_sequential_lin_and_assign(obj, map, pattern, ln_prefact, struct('display', 1, 'maxit', 30), @(x) assign_perm(x, pattern{1}, [1, 0, 1, 1]), 1);

        pattern = {[a, 0, c, b]};

        obj = solve_sequential_lin_and_assign(obj, map, pattern, ln_prefact, struct('display', 1, 'maxit', 30), @(x) assign_perm(x, pattern{1}, [1, 0, 1, 1]), 1);

        %alpha_pattern_perm = {rot_90};
        %[obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);

        err1 = calculate_error(obj, map, [], 1);

        if obj.testing == 1
            if err1 > obj.copts.err_tol
                disp(err1);
            end
        end
    end
    %obj = cell2matrix(obj); max(reshape( obj.PEPO_matrix ,[],1))

    %     if obj.testing == 1
    %         calculate_error(obj, [
    %                         1, 1, 1;
    %                         1, 1, 1; ], [], 1)
    %         calculate_error(obj, [
    %                         1, 1, 0;
    %                         1, 1, 1;
    %                         0, 1, 1; ], [], 1)
    %         %         calculate_error(obj, [
    %         %                         1, 1, 1;
    %         %                         1, 1, 1;
    %         %                         0, 1, 1; ], [], 1)
    %     end

end
