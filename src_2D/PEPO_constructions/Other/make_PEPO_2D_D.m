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

    L = 2;
    max_bond_dim = 16;

    for n = 1:L
        sz = min(d^(2 * n), max_bond_dim);

        obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, sz];
        obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, sz];
    end
    %
    %     a = numel(obj.virtual_level_sizes_horiz);
    %     b = a+1;
    %     %b = a + 1;
    %     c = b + 1;
    %     %e = c + 1;
    %     l=c+1;
    %
    %     obj.current_max_index = l+1;
    %     obj.max_index = l+1;
    %
    %     dd = [8,8,6,d^2,d^2];

    %
    l = numel(obj.virtual_level_sizes_horiz);
    a = l + 2;
    b = a;

    dd = [d^2, d^2, 10];

    %dd = [8, 8, 6];
    %dd = [10, 10, 6, 6];

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, dd];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, dd];

    %% function definitions for adding linear and loop blocks

    function [obj, ln_prefact, err] = add_lin(obj, pattern, ln_prefact, nl, cons)
        if nargin < 4
            nl = 0;
        end

        if nargin < 5
            cons = 0;
        end

        [map1, ~] = create_map(make_cross(pattern));

        if cons == 1
            c2 = [pattern ~= 0] * 2;
            pat2 = pattern;
            pat2(pattern ~= 0) = pat2(pattern ~= 0) + l -1;
            pat3 = pat2 + [pattern ~= 0] .* [0, 0, 1, 1];
        else
            pat2 = pattern;
            pat3 = pattern;

            c2 = [0, 0, 0, 0];
        end

        if nl == 0
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, {pat3}, ln_prefact);
        else
            [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map1}, {pat3}, ln_prefact, nl_opts, rot_180);
        end
        obj = assign_perm(obj, pat2, c2);

        err = calculate_error(obj, map1, [], 1);

        if obj.testing == 1
            if err > obj.err_tol
                disp(err);
            end
        end
    end

    function [obj, ln_prefact] = add_loop(obj, arr, pat, ln_prefact)
        [map1, ~] = make_cross_loop(arr, [b, a]);
        map1 = create_map(map1, struct);

        pat2 = [pat, b, a];

        %if nl == 0
        [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, {pat2}, ln_prefact);
        %else
        %    [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map1}, {pat}, ln_prefact, nl_opts, rot_180);
        %end

        obj = assign_perm(obj, pat2, [2 * (pat == l), 1, 1]);

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

    %     obj.PEPO_cell{1, 1, l+1, 1 } = reshape(  eye(d^2) , [d, d, 1,1,d^2,  1]);
    %     obj = assign_perm(obj, [0, 0, l, 0], [0, 0, -2, 0]);
    %
    %
    %     [obj, ln_prefact, ~] = add_lin(obj, [1, 1, 1, 0], ln_prefact,0,1);
    %     [obj, ln_prefact, ~] = add_lin(obj, [1, 1, 1, 1], ln_prefact,0,1);

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
                [obj, ln_prefact, ~] = add_lin(obj, [3, 3, 3, 3], ln_prefact);
            otherwise
                error('not impl')
        end
    end

    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% LOOPS %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %     if obj.testing == 1
    %         fprintf('started loops')
    %     end
    %
    %simple loop
    [map, ~] = create_map([1, 2; 3, 4], obj.numopts);

    alpha_pattern_perm = {rot_90, rot_90, rot_90, rot_90};

    pattern = {[a, a, 0, 0]};

    [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);

    %
    obj.PEPO_cell{1, 1, l + 1, 1} = reshape(eye(d^2) / exp(obj.nf), [d, d, 1, 1, d^2, 1]);

    map = create_map([1, 1, 1;
                0, 1, 1], struct);

    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, {[l, 0, a, a]}, ln_prefact);
    obj = rescale_PEPO_pattern(obj, {[0, 0, l, 0], [l, 0, a, a]});

    obj = assign_perm(obj, [l, 0, a, a], [2, 0, 1, 1]);
    obj = assign_perm(obj, [0, 0, l, 0], [0, 0, -2, 0]);

    %
    [obj, ln_prefact] = add_loop(obj, [1, 0], [l, 0], ln_prefact);
    [obj, ln_prefact] = add_loop(obj, [1, 1], [l, l], ln_prefact);

    %
    %     if obj.testing == 1
    %         fprintf('double loops')
    %     end
    %

    %
    %     if obj.testing == 1
    %         fprintf('started loops')
    %     end
    %     %simple loop
    %     [map, ~] = create_map([1, 2; 3, 4], obj.numopts);
    %
    %     alpha_pattern_perm = {rot_90, rot_90, rot_90};
    %     pattern = {[c, c, 0, 0]};
    %
    %     [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);
    %
    %     %
    %     pattern = {[b, 0, 0, c], ...
    %                 [0, a, c, 0],...
    %                 [0, 0, l, 0]};
    %
    %     obj = fill_rand(obj, pattern);
    %
    %     [map, ~] = create_map([
    %                         1, 1, 1;
    %                         0, 1, 1; ], struct);
    %     pattern2 = {[l, 0, b, a]};
    %     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern2, ln_prefact);
    %
    %     obj = rescale_PEPO_pattern(obj, [pattern, pattern2]);
    %
    %     obj = assign_perm(obj, [0, 0, l, 0], [0, 0, -2, 0]);
    %     obj = assign_perm(obj, [l, 0, b, a], [2, 0, 1, 1]);
    %     obj = assign_perm(obj, [b, 0, 0, c], [1, 0, 0, 1]);
    %     obj = assign_perm(obj, [0, a, c, 0], [0, 1, 1, 0]);
    %
    %    %calculate 2nd time
    %
    %     [obj, ln_prefact] = add_loop(obj, [1, 0],[l, 0], ln_prefact);
    %     [obj, ln_prefact] = add_loop(obj, [1, 1],[l, l], ln_prefact);
    %

end

function [m, p] = rotate(m, p, k)
    m = rot90(m, k);
    for l = 1:numel(p)
        p{l} = circshift(p{l}, -k);
    end
end
