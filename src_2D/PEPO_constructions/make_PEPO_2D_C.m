function [obj, error_code] = make_PEPO_2D_C(obj)

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
        nl_opts = struct('Gradient', true, 'Display', 'iter', 'maxit', 20);
    else
        nl_opts = struct('Gradient', true, 'Display', 'none', 'maxit', 20);
    end

    obj.current_max_index = 1;
    obj.virtual_level_sizes_horiz = [1];
    obj.virtual_level_sizes_vert = [1];

    obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
    obj.bounds = [1];
    obj.boundary_vect(obj.bounds) = 1;

    obj.complex = false;

    L = 2;
    max_bond_dim = 16;

    for n = 1:L
        sz = min(d^(2 * n), max_bond_dim);

        obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, sz];
        obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, sz];
    end

    a = numel(obj.virtual_level_sizes_horiz);
    b = a;
    %b = a + 1;
    %c = b + 1;
    %e = c + 1;

    obj.current_max_index = a;
    obj.max_index = a;

    dd = [6];

    %dd = [8, 8, 6];
    %dd = [10, 10, 6, 6];

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
        [map1, pat] = make_cross_loop(arr, [b, a]);
        map1 = create_map(map1, struct);
        %if nl == 0
        [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, {pat}, ln_prefact, struct);
        %else
        %    [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map1}, {pat}, ln_prefact, nl_opts, rot_180);
        %end

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
                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 3, 2], ln_prefact);
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

    alpha_pattern_perm = {rot_90, rot_90, rot_90, rot_90};
    %     pattern = {...
    %                 [0, 0, b, a] ...
    %                 [b, 0, 0, c], ...
    %                 [e, c, 0, 0], ...
    %                 [0, a, e, 0]};

    pattern = {[a, a, 0, 0]};

    [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);

    %     pattern = {[a,a,0,0]};
    %
    %     [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);
    %
    %     tt = obj.PEPO_cell{a+1,a+1,1,1}/ (4^(1/4));
    %
    %     obj.PEPO_cell{b+1,a+1,1,1} = tt;
    %     obj.PEPO_cell{c+1,b+1,1,1} = tt;
    %     obj.PEPO_cell{e+1,c+1,1,1} = tt;
    %     obj.PEPO_cell{a+1,e+1,1,1} = tt;
    %
    %     obj.PEPO_cell{a+1,a+1,1,1}=[];
    %     obj.PEPO_cell{1,a+1,a+1,1}=[];
    %     obj.PEPO_cell{1,1,a+1,a+1}=[];
    %     obj.PEPO_cell{a+1,1,1,a+1}=[];
    %
    %     obj = assign_perm(obj, [b,a,0,0], [1, 1, 0, 0]);
    %     obj = assign_perm(obj, [c,b,0,0], [1, 1, 0, 0]);
    %     obj = assign_perm(obj, [e,c,0,0], [1, 1, 0, 0]);
    %     obj = assign_perm(obj, [a,e,0,0], [1, 1, 0, 0]);

    %loops with 1 or more legs
    %     for n = 1:L
    %         switch n
    %             case 1
    %                 [obj, ln_prefact] = add_loop(obj, [1, 0], ln_prefact);
    %                 [obj, ln_prefact] = add_loop(obj, [1, 1], ln_prefact);
    %             case 2
    %                 [obj, ln_prefact] = add_loop(obj, [2, 0], ln_prefact);
    %                 [obj, ln_prefact] = add_loop(obj, [2, 1], ln_prefact);
    %                 [obj, ln_prefact] = add_loop(obj, [2, 2], ln_prefact);
    %             case 3
    %                 [obj, ln_prefact] = add_loop(obj, [3, 0], ln_prefact);
    %                 [obj, ln_prefact] = add_loop(obj, [3, 1], ln_prefact);
    %                 [obj, ln_prefact] = add_loop(obj, [3, 2], ln_prefact);
    %                 [obj, ln_prefact] = add_loop(obj, [3, 3], ln_prefact);
    %             otherwise
    %                 error('not impl')
    %         end
    %     end
    %
    %     if obj.testing == 1
    %         fprintf('double loops')
    %     end

    %offset loop
    %     [map, ~] = create_map([
    %                         1, 1, 0;
    %                         1, 1, 1;
    %                         0, 1, 1; ], struct);
    %
    %
    %     pattern = {[a, a, a, a]};
    %     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %      obj = assign_perm(obj, pattern{1}, [0, 0, 0, 0]);

    %double loop
    %
    %           alpha_pattern_perm = {rot_90};
    %           [map, ~] = create_map([
    %                               1, 1, 1;
    %                               1, 1, 1; ], struct);
    %           pattern = {[a, 0, a, a]};
    %
    %           %obj = solve_sequential_lin_and_assign(obj, map, pattern, ln_prefact, struct('display', 1), @(x) assign_perm(x, [a,0,a,a], [1, 0, 1, 1]), 0.5)
    %
    %           [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);

end

function [m, p] = rotate(m, p, k)
    m = rot90(m, k);
    for l = 1:numel(p)
        p{l} = circshift(p{l}, -k);
    end
end
