function [obj, error_code] = make_PEPO_2D_A(obj)
    d = obj.dim;
    ln_prefact = obj.nf;

    error_code = 0;

    rot_180 = {{[3, 4, 1, 2]}};
    if obj.testing == 1
        nl_opts = struct('Gradient', true, 'Display', 'iter', 'maxit', 100);
    else
        nl_opts = struct('Gradient', true, 'Display', 'none', 'maxit', 100);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Block with 1/2 legs %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%% LEVEL 1 %%%%%%%%%%%%%%
    obj.current_max_index = 1;
    obj.virtual_level_sizes_horiz = [1];
    obj.virtual_level_sizes_vert = [1];

    obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
    obj.bounds = [1];
    obj.boundary_vect(obj.bounds) = 1;

    obj.PEPO_cell{1, 1, 1, 1} = reshape((expm(obj.H_1_tensor)) / exp(obj.nf), [d, d, 1, 1, 1, 1]);

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

        obj = cell2matrix(obj);
        err = calculate_error(obj, map1, [], 1);

        if obj.testing == 1
            if err > obj.err_tol
                disp(err);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Blocks 1D chain %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    obj.complex = true;

    L = 2;

    for l = 1:L
        sz = min(d^(2 * l), Inf);

        obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, sz];
        obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, sz];
    end

    for l = 1:L

        %%% n-1--|--n--|--n1- block

        if l ~= 1
            err01 = calculate_error(obj, 1:(2 * l + 2), obj.cycleopts, 1);
        else
            err01 = 1;
        end

        [obj, ln_prefact, ~] = add_lin(obj, [l - 1, 0, l, 0], ln_prefact, 1);
        err02 = calculate_error(obj, 1:(2 * l + 2), obj.cycleopts, 1);

        if err02 > err01
            warning('not converging')
            error_code = 1;

            obj.PEPO_cell{l + 1, l + 1, 1, 1} = [];
            obj = assign_perm(obj, [l, l, 0, 0]);
        end

        %%% n--|--n block

        [obj, ln_prefact, ~] = add_lin(obj, [l, l, 0, 0], ln_prefact);
        err03 = calculate_error(obj, 1:(2 * l + 2), obj.cycleopts, 1);

        if err03 -1e-14 > err02
            warning('not converging')
            error_code = 1;
            obj.PEPO_cell{l + 1, l + 1, 1, 1} = [];
            obj = assign_perm(obj, [l, l, 0, 0]);

        end

        % blocks with 3 or more legs

        switch l
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

        %make loop correction at current stage
        %         switch l
        %             case 1
        %                 [map, ~] = create_map([1, 2; 3, 4], obj.numopts);
        %
        %                 pattern = {[0, 0, a, 1],...
        %                              [a,0,0,1]};
        %
        %                 rpat={[2, 3, 4, 1], ...
        %                 [3, 4, 1, 2], ...
        %                 [4, 1, 2, 3]};
        %
        %                 [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, {rpat,rpat});
        %
        %                 %
        %                 [map, ~] = create_map([1, 1, 1;
        %                                        0, 1, 1], struct);
        %                 pattern = {[1, 0, a, 1]};
        %
        %                 %[obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
        %
        %                 %[obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts);
        %                 obj = assign_perm(obj, [1, 0, a, 1], [0,0,1,1] );
        %                 %
        % %                 [map, ~] = create_map([1, 1, 1;
        % %                                        1, 1, 0], struct);
        % %                 pattern = {[a, 0, 1, 1]};
        % %                 [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts);
        % %                 %
        %                 [map, ~] = create_map([0, 1, 0;
        %                                        1, 1, 1;
        %                                        0, 1, 1], struct);
        %                 pattern = {[1, 1, a, 1]};
        %                 [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts);
        %                 obj = assign_perm(obj, [1, 1, a, 1] , [0,0,1,1] );
        %
        %             case 2
        %
        %             case 3
        %
        %         end

        %         [map, ~] = create_map([1, 2; 3, 4], obj.numopts);
        %
        %         rpat = {...
        %                 [2, 3, 4, 1], ...
        %                 [3, 4, 1, 2], ...
        %                 [4, 1, 2, 3]};
        %         if l ==1
        %             pattern = {[0, 0, l, a],...
        %                         [b,l,0,0],...
        %                         [0,a,b,0]};
        %
        %             alpha_pattern_perm = {rpat,rpat,rpat};
        %         else
        %             pattern = {[0, 0, l, a],...
        %                         [b,l,0,0]};
        %             alpha_pattern_perm = {rpat,rpat};
        %         end
        %
        %         [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts,alpha_pattern_perm);
    end

    %%
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% LOOPS %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    a = numel(obj.virtual_level_sizes_horiz);
    b = a + 1;
    c = b + 1;

    obj.current_max_index = c;
    obj.max_index = c;

    dd = [d^3, d^3, d];

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, dd];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, dd];

    %d0 = numel(obj.virtual_level_sizes_horiz);
    %d1 = d0+1;
    %a = d0+2;

    %     a = numel(obj.virtual_level_sizes_horiz);
    %     b=a+1;
    %
    %     obj.current_max_index = b;
    %     %obj.PEPO_cell{zero_level+1,zero_level+1,zero_level+1,zero_level+1}= obj.PEPO_cell{1,1,1,1};
    %
    %     alpha_dim = 6;
    %
    %     obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, alpha_dim, alpha_dim];
    %     obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, alpha_dim , alpha_dim];
    %
    nl_opts = struct('Gradient', true, 'Display', 'iter');

    %      [map, ~] = create_map([1, 2; 3, 4], obj.numopts);
    %      pattern = {[a, a, 0, 0]};
    %
    rpat = {[2, 3, 4, 1], ...
            [3, 4, 1, 2], ...
            [4, 1, 2, 3]};

    alpha_pattern_perm = {rpat, rpat, rpat, rpat};
    %
    %      [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);

    [map, ~] = create_map([1, 2; 3, 4], obj.numopts);
    pattern = {[0, 0, b, a], ...
                [b, 0, 0, c], ...
                [c, c, 0, 0], ...
                [0, a, c, 0]};

    [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);

    calculate_error(obj, map, [])

    %single extension
    [map, ~] = create_map([
                        0, 1, 0, 0;
                        0, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 0, 0, 0]);
    pattern = {[0, 1, b, a]};

    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1}, [0, 0, 1, 1]);

    %double extension, can be solve linearly
    [map, ~] = create_map([
                        0, 1, 0, 0;
                        1, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 0, 0, 0]);

    pattern = {[1, 1, b, a]};

    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1}, [0, 0, 1, 1]);

    %long extension
    [map, ~] = create_map([
                        0, 0, 1, 0;
                        0, 0, 1, 0;
                        0, 0, 1, 1;
                        0, 0, 1, 1]);
    pattern = {[0, 2, b, a]};

    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    obj = assign_perm(obj, pattern{1}, [0, 0, 1, 1]);

    [map, ~] = create_map([
                        0, 0, 1, 0;
                        0, 0, 1, 0;
                        0, 1, 1, 1;
                        0, 0, 1, 1]);
    pattern = {[1, 2, b, a]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1}, [0, 0, 1, 1]);
    %
    [map, ~] = create_map([
                        0, 0, 1, 0;
                        0, 0, 1, 0;
                        1, 1, 1, 1;
                        0, 0, 1, 1;
                        ]);
    pattern = {[2, 2, b, a]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1}, [0, 0, 1, 1]);
    %
    if obj.testing == 1
        obj = cell2matrix(obj);
        calculate_error(obj, [
                        0, 0, 1, 0;
                        0, 0, 1, 0;
                        0, 0, 1, 1;
                        0, 0, 1, 1;
                        ], struct, 1)
    end

    %double 2 options

    %      for k = 0:3
    %
    %         [map, pattern] = rotate([
    %                             0, 0, 1, 0;
    %                             1, 1, 1, 0;
    %                             0, 1, 1, 0;
    %                             ], {
    %         [1, 0, b, a], [b, 1, 0, a]
    %         }, k);
    %         [map, ~] = create_map(map);
    %         obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnopts);
    %
    %         [map, pattern] = rotate([
    %                             0, 0, 0, 0;
    %                             1, 1, 1, 1;
    %                             0, 1, 1, 0;
    %                             ], {
    %         [b, 0, 1, a]
    %         }, k);
    %         [map, ~] = create_map(map);
    %         obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnopts);
    %
    %         [map, pattern] = rotate([
    %                             0, 0, 1, 0;
    %                             1, 1, 1, 1;
    %                             0, 1, 1, 0;
    %                             ], {
    %         [b, 1, 1, a]
    %         }, k);
    %         [map, ~] = create_map(map);
    %         obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnopts);
    %
    %         [map, pattern] = rotate([
    %                             0, 1, 0, 0;
    %                             1, 1, 1, 1;
    %                             0, 1, 1, 0;
    %                             ], {
    %         [1, 1, b, a]
    %         }, k);
    %         [map, ~] = create_map(map);
    %         obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnopts);
    %
    %         [map, pattern] = rotate([
    %                             0, 0, 1, 0;
    %                             1, 1, 1, 1;
    %                             0, 1, 1, 0;
    %                             ], {
    %         [b, 1, 1, a]
    %         }, k);
    %         [map, ~] = create_map(map);
    %         obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnopts);
    %
    %     end

    %
    %     if obj.testing == 1
    %         obj = cell2matrix(obj);
    %         calculate_error(obj, [
    %                         0, 0, 1, 0;
    %                         1, 1, 1, 0;
    %                         0, 1, 1, 0;
    %                         0, 0, 0, 0], struct, 1)
    %     end

    %fprintf(".");
end

% function error_code = check_lin_error(obj)
%     err = calculate_error(obj, 1:6, struct("numbered", true, "h_cyclic", 1));
%     error_code = err > 1 + 1e-3;
% end

function [m, p] = rotate(m, p, k)
    m = rot90(m, k);
    for l = 1:numel(p)
        p{l} = circshift(p{l}, -k);
    end
end

function test_rot_90(obj)
    obj = cell2matrix(obj);
    svds(reshape(obj.PEPO_matrix - permute(obj.PEPO_matrix, [1, 2, 5, 6, 3, 4]), [], 1))

end
