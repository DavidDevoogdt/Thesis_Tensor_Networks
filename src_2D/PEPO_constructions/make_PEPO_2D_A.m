function [obj, error_code] = make_PEPO_2D_A(obj)

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

    obj.current_max_index = 1;
    obj.virtual_level_sizes_horiz = [1];
    obj.virtual_level_sizes_vert = [1];

    obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
    obj.bounds = [1];
    obj.boundary_vect(obj.bounds) = 1;

    obj.complex = false;

    L = obj.L;
    max_bond_dim = obj.max_bond_dim;

    for n = 1:L
        sz = min(d^(2 * n), max_bond_dim);

        obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, sz];
        obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, sz];
    end

    a = numel(obj.virtual_level_sizes_horiz);
    b = a + 1;
    c = b + 1;

    obj.current_max_index = c;
    obj.max_index = c;

    %dd = [8, 8, 6];
    dd = [6,10, 10];

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, dd];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, dd];

    %% function definitions for adding linear and loop blocks

    function [obj, ln_prefact, err] = add_lin(obj, pattern, ln_prefact, nl)
        if nargin < 4
            nl = 0;
        end

        [map1, ~] = create_map(make_cross(pattern));

        if nl == 0
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, {pattern}, ln_prefact, struct);
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

    function [obj, ln_prefact] = add_loop(obj, arr, ln_prefact)
        pat = [arr, c, b];

        map1 = make_cross_loop(pat, [0, 0, 1, 1]);
        map1 = create_map(map1, struct);

        [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, {pat}, ln_prefact, struct);

        obj = assign_perm(obj, pat, [0, 0, 1, 1]);

        err = calculate_error(obj, map1, [], 1);

        if obj.testing == 1
            if err > obj.err_tol
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

        [obj, ln_prefact, ~] = add_lin(obj, [n - 1, 0, n, 0], ln_prefact, 1);
        err02 = calculate_error(obj, 1:(2 * n + 2), obj.cycleopts, 1);

        if err02 > err01
            warning('not converging')
            error_code = 1;

            obj.PEPO_cell{n + 1, n + 1, 1, 1} = [];
            obj = assign_perm(obj, [n, n, 0, 0]);
        end

        %%% n--|--n block

        [obj, ln_prefact, ~] = add_lin(obj, [n, n, 0, 0], ln_prefact);
        err03 = calculate_error(obj, 1:(2 * n + 2), obj.cycleopts, 1);

        if err03 > err02
            warning('not converging')
            error_code = 1;
            obj.PEPO_cell{n + 1, n + 1, 1, 1} = [];
            obj = assign_perm(obj, [n, n, 0, 0]);
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
                [obj, ln_prefact, ~] = add_lin(obj, [1, 1, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [1, 1, 1, 1], ln_prefact);

            case 2
                [obj, ln_prefact, ~] = add_lin(obj, [2, 1, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [2, 1, 1, 1], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [2, 2, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [2, 2, 1, 1], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [2, 2, 2, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [2, 2, 2, 1], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [2, 2, 2, 2], ln_prefact);
            case 3
                [obj, ln_prefact, ~] = add_lin(obj, [3, 1, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 1, 1, 1], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [3, 2, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 2, 1, 1], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [3, 2, 2, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 2, 2, 1], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 2, 2, 2], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 1, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 1, 1], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 2, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 2, 1], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 2, 2], ln_prefact);

                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 3, 0], ln_prefact);
                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 3, 1], ln_prefact);

                %works but slow
                %[obj,ln_prefact,~] = add_lin(obj, [3,3,3,2] ,ln_prefact);
                %[obj,ln_prefact,~] = add_lin(obj, [3,3,3,3] ,ln_prefact);
            otherwise
                error('not impl')
        end
    end

    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% LOOPS %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if obj.testing == 1
        fprintf('started loops')
    end
    %simple loop
    [map, ~] = create_map([1, 2; 3, 4], obj.numopts);

    alpha_pattern_perm = {rot_90, rot_90, rot_90};
    pattern = {[a, a, 0, 0]};

    [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);

    loop_extension = 1;
    
    if loop_extension==1
    
        
        
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
        [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern2, ln_prefact, struct);

        obj = rescale_PEPO_pattern(obj, [pattern, pattern2]);

        obj = assign_perm(obj, [1, 0, c, b], [0, 0, 1, 1]);
        obj = assign_perm(obj, [c, 0, 0, a], [1, 0, 0, 1]);
        obj = assign_perm(obj, [0, b, a, 0], [0, 1, 1, 0]);

        err1 = calculate_error(obj, map, [], 1);

        if obj.testing == 1
            if err1 > obj.err_tol
                disp(err1);
            end
        end
        
        %loops with 1 or more legs
        for n = 1:L
            switch n
                case 1
                    [obj, ln_prefact] = add_loop(obj, [1, 0], ln_prefact);
                    [obj, ln_prefact] = add_loop(obj, [1, 1], ln_prefact);
                case 2
                    [obj, ln_prefact] = add_loop(obj, [2, 0], ln_prefact);
                    [obj, ln_prefact] = add_loop(obj, [2, 1], ln_prefact);
                    [obj, ln_prefact] = add_loop(obj, [2, 2], ln_prefact);
                case 3
                    [obj, ln_prefact] = add_loop(obj, [3, 0], ln_prefact);
                    [obj, ln_prefact] = add_loop(obj, [3, 1], ln_prefact);
                    [obj, ln_prefact] = add_loop(obj, [3, 2], ln_prefact);
                    [obj, ln_prefact] = add_loop(obj, [3, 3], ln_prefact);
                otherwise
                    error('not impl')
            end
        end
    
    end

    

    %     [map, ~] = create_map([
    %                         1, 1, 0;
    %                         1, 1, 1;
    %                         0, 1, 1; ], struct);
    %     %pattern = {[b, a, b, a]};
    %     pattern = {[c, a, b, c]};
    %     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    %     obj = scale_PEPO_pattern(obj, pattern, 0.5 );
    %     obj = assign_perm(obj, pattern{1}, [1, 1, 1, 1]);
    %
    %     err1 = calculate_error(obj, map, [], 1);
    %
    %     if obj.testing == 1
    %         if err1 > obj.err_tol
    %             disp(err1);
    %         end
    %     end

    %double loop

    %    alpha_pattern_perm = {rot_90};
    %    [map, ~] = create_map([
    %                        1, 1, 1;
    %                        1, 1, 1; ], struct);
    %    pattern = {[c, 0, c, c]};
    %
    %    %obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, struct('display', 1), @(x) assign_perm(x, [a,0,a,a], [1, 0, 1, 1]), 0.1)
    %
    %    [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);
    %
    %     err1 = calculate_error(obj, map, [], 1);
    %
    %     if obj.testing == 1
    %         if err1 > obj.err_tol
    %             disp(err1);
    %         end
    %     end

    %obj = cell2matrix(obj); max(reshape( obj.PEPO_matrix ,[],1))

    if obj.testing == 1
        calculate_error(obj, [
                        1, 1, 1;
                        0, 1, 1; ], [], 1)
        calculate_error(obj, [
                        1, 1, 0;
                        1, 1, 1;
                        0, 1, 1; ], [], 1)
        %         calculate_error(obj, [
        %                         1, 1, 1;
        %                         1, 1, 1;
        %                         0, 1, 1; ], [], 1)
    end

end
