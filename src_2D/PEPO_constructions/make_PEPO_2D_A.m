function [obj, error_code ] = make_PEPO_2D_A(obj)
    d = obj.dim;
    ln_prefact = obj.nf;

    error_code = 0;

    rot_180 = {{[3, 4, 1, 2]}};

    %nl_opts = struct('Gradient', true, 'Display', 'iter-detailed','Algoritm', "trust-region");
    nl_opts = struct('Gradient', true, 'Display', 'None');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Block with 1/2 legs %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%% LEVEL 1 %%%%%%%%%%%%%%
    obj.current_max_index = 1;
    obj.virtual_level_sizes_horiz = [1, d^2];
    obj.virtual_level_sizes_vert = [1, d^2];

    obj.boundary_vect = zeros(1, size(obj.PEPO_cell, 1));
    obj.bounds = [1];
    obj.boundary_vect(obj.bounds) = 1;

    obj.PEPO_cell{1, 1, 1, 1} = reshape((expm(obj.H_1_tensor)) / exp(obj.nf), [d, d, 1, 1, 1, 1]);

    %obj.PEPO_cell{1,1,1,1} = reshape(eye(d) / exp(obj.nf), [d, d, 1, 1, 1, 1]);
    %obj.PEPO_cell{one_level+1,one_level+1,one_level+1,one_level+1}  =  reshape( (expm(obj.H_1_tensor)- eye(d) ) / exp(obj.nf), [d, d, 1, 1, 1, 1]);

    
    %0--|--1--|--0 and all other veriants
    for i=1:2
        n = 2;
        [map, ~] = create_map(1:n, obj.numopts);
        pattern = {[1, 0, 0, 0]};
        [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, rot_180);

        obj = assign_perm(obj, pattern{1});
        
        err = calculate_error(obj, 1:n, obj.numopts);
        if err > obj.err_tol
            if i ==1
                fprintf("switching to complex tensors")
                obj.complex = true;
            else
                if obj.testing == 1
                    disp(err);
                end 
            end
        else
           break; 
        end
    
    end
    

    % 0--|--1--|1--|--0 and all other veriants
    %#6
    n = 3;
    [map, ~] = create_map(1:n, obj.numopts);
    pattern = {[1, 0, 1, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    %eq_perm = {[1, 2, 4, 3]};
    obj = assign_perm(obj, [1, 0, 1, 0]);

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, [1 2 3], obj.numopts);
        if err > obj.err_tol
            disp(err);
        end
        
        %         calculate_error(obj, [1 2; 0 3], obj.numopts)
        %         calculate_error(obj, [1 0; 2 3], obj.numopts)
        %         calculate_error(obj, [1; 2; 3; ], obj.numopts)
        %
        %         calculate_error(obj, [1, 2; 3, 0], obj.numopts)
        %         calculate_error(obj, [0, 2; 3, 1], obj.numopts)
    end

    %%
    %%%%%%%%%%%%%% LEVEL 2 %%%%%%%%%%%%%%
    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^4];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^4];
    obj.current_max_index = 2;
    %%
    %%%%--1--|--2--|--1 and variants%%%%
    %# 4*3
    [map, ~] = create_map(1:4, obj.numopts);

    pattern = {[1, 0, 2, 0]};
    [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, rot_180);

    obj = assign_perm(obj, [1, 0, 2, 0]);

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
        %         calculate_error(obj, [0, 0, 4; 1, 2, 3], obj.numopts)
        %         calculate_error(obj, [1, 2, 4; 3, 0, 0], obj.numopts)
        %         calculate_error(obj, [1, 2; 3, 0; 4, 0], obj.numopts)
        %         calculate_error(obj, [0, 2; 0, 1; 3, 4], obj.numopts)
        %         calculate_error(obj, (1:4)', obj.numopts)
    end

    %%
    %%%%--2--|--2-- and variants%%%%
    %# 4*3/2= 6

    [map, ~] = create_map(1:5, obj.numopts);
    pattern = {[2, 0, 2, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    obj = assign_perm(obj, [2, 0, 2, 0]);

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        %calculate_error(obj, [1 2 3; 0 0 4; 0 0 5], obj.numopts)
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
        %         calculate_error(obj, [1 2 3 0; 0 0 4 5], obj.numopts)
        %         calculate_error(obj, [2, 3, 4, 1; 0, 0, 0, 5], obj.numopts)
        %         calculate_error(obj, [2, 3, 4, 1; 5, 0, 0, 0], obj.numopts)
        %
        %         calculate_error(obj, [0, 3, 4, 5; 1, 2, 0, 0], obj.numopts)
        %         calculate_error(obj, [0, 0, 4, 5; 1, 2, 3, 0], obj.numopts)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Block with 3/4 legs %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    obj.current_max_index = obj.max_index;

    %%%%%%%%%%%%%% LEVEL 1 %%%%%%%%%%%%%%
    %% 1 block with 3 external 1 legs

    [map, ~] = create_map([
                        0, 1, 0;
                        1, 1, 1]);
    pattern = {[1, 1, 1, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    obj = assign_perm(obj, pattern{1});

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
        %         calculate_error(obj, [3, 2, 4; 0, 1, 0], obj.numopts)
        %         calculate_error(obj, [1, 0; 2, 3; 4, 0], obj.numopts)
        %         calculate_error(obj, [0, 1; 2, 3; 0, 4], obj.numopts)
    end
    %% 1 block with 4 external legs
    [map, ~] = create_map([
                        0, 1, 0;
                        1, 1, 1;
                        0, 1, 0]);
    pattern = {[1, 1, 1, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1});

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
    end

    %%%%%%%%%%%%%% LEVEL 2 %%%%%%%%%%%%%%

    %% 1 2 a and 2 1 levels

    [map, ~] = create_map([
                        1, 1, 1, 1;
                        0, 1, 0, 0]);
    pattern = {[1, 0, 2, 1]};

    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    obj = assign_perm(obj, pattern{1});

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
        %calculate_error(obj, [1, 0; 2, 0; 4, 3; 5, 0], obj.numopts)
    end

    %% 1 2 a and 3 1 levels

    [map, ~] = create_map([
                        0, 1, 0, 0;
                        1, 1, 1, 1;
                        0, 1, 0, 0]);
    pattern = {[1, 1, 2, 1]};
    %eq_perm = {};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1});

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
    end

    %% 2 2 levels + 1 a

    [map, ~] = create_map([
                        1, 1, 1, 1, 1;
                        0, 0, 1, 0, 0]);
    pattern = {[2, 0, 2, 1]};
    %eq_perm = {[2, 4, 1, 3], [4, 2, 1, 3]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1});

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
        %calculate_error(obj, [0, 1; 0, 2; 3, 4; 0, 5; 0, 6], obj.numopts)

    end

    %% 2 2 levels + 2 1 a

    [map, ~] = create_map([
                        0, 0, 1, 0, 0;
                        1, 1, 1, 1, 1;
                        0, 0, 1, 0, 0]);
    pattern = {[2, 1, 2, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1});

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
    end

    %%

    [map, ~] = create_map([
                        1, 1, 1, 1, 1;
                        0, 0, 1, 0, 0;
                        0, 0, 1, 0, 0; ]);
    pattern = {[2, 0, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1});

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
        %calculate_error(obj, [0, 0, 7, 0, 0; 0, 0, 3, 0, 0; 1, 2, 6, 4, 5], obj.numopts)

    end

    % 3 2 levels + 1 1 level

    [map, ~] = create_map([
                        0, 0, 1, 0, 0;
                        1, 1, 1, 1, 1;
                        0, 0, 1, 0, 0;
                        0, 0, 1, 0, 0; ]);
    pattern = {[2, 1, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1});

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
        %calculate_error(obj, [0, 0, 7, 0, 0; 0, 0, 3, 0, 0; 1, 2, 6, 4, 5], obj.numopts)

    end

    %%
    %4 2
    [map, ~] = create_map([
                        0, 0, 1, 0, 0;
                        0, 0, 1, 0, 0;
                        1, 1, 1, 1, 1;
                        0, 0, 1, 0, 0;
                        0, 0, 1, 0, 0]);
    pattern = {[2, 2, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    obj = assign_perm(obj, pattern{1});

    error_code = check_lin_error(obj);
    if error_code ~= 0
       return ; 
    end
    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
    end
    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% LOOPS %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %x = 2;
    x = 2;
    a = 3;
    b = a + 1;

    %obj.PEPO_cell{zero_level+1,zero_level+1,zero_level+1,zero_level+1}= obj.PEPO_cell{1,1,1,1};

    alpha_dim = 8;
    %beta_dim = 20;

    %obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^2, alpha_dim];
    %obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^2, alpha_dim];

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, alpha_dim];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, alpha_dim];

    obj.current_max_index = numel(obj.virtual_level_sizes_horiz);
    obj.max_index = obj.current_max_index;

    %simple loop

    [map, ~] = create_map([1, 2; 3, 4], obj.numopts);

    pattern = {[a, a, 0, 0]};
    alpha_pattern_perm = {{...
                            [2, 3, 4, 1], ...
                            [3, 4, 1, 2], ...
                            [4, 1, 2, 3]}};

    %nl_opts = struct('Gradient', true, 'Display', 'iter-detailed');

    lnopts = struct('Display', 0, 'maxit', 1);
    [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);
    %obj = assign_perm(obj, alpha_pattern{1},a);

    
    if obj.testing == 1
        err = calculate_error(obj, map);
        if err > obj.err_tol
            disp(err);
        end
    end

    obj = cell2matrix(obj);
    err = calculate_error(obj, [
                             1, 1, 1;
                             1, 1, 1;], struct, 1);
                         
    if err > 1
       error_code = 1;
    end
    
    
    
%     %nl_opts = struct('Gradient', true, 'Display', 'iter-detailed', 'Algoritm', "trust-region");
%      nl_opts = struct('Gradient', true, 'Display', 'iter-detailed');
% 
% 
%         [map, ~] = create_map([1, 2; 3, 4; 5, 6], obj.numopts);
%         pattern = {[a, a, 0, a]};
%         pattern_perm = {{...
%                                 [3, 2, 1, 4]}};
%         [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, pattern_perm);
% 
%         if obj.testing == 1
%             calculate_error(obj, [1, 2; 3, 4; 5, 6], obj.numopts)
%          end

    % %% double loops

    %         %single extension
%             [map, ~] = create_map([
%                                 0, 1, 0, 0;
%                                 0, 1, 1, 0;
%                                 0, 1, 1, 0;
%                                 0, 0, 0, 0]);
%             pattern = {[0, 1, a, a]};
%             %[obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
%             [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts);
%             obj = assign_perm(obj, pattern{1}, [0, 0, 1, 1]);
    %
    %         %double extension
%             [map, ~] = create_map([
%                                 0, 1, 0, 0;
%                                 1, 1, 1, 0;
%                                 0, 1, 1, 0;
%                                 0, 0, 0, 0]);
%             pattern = {[ 1,  1, a, a]};
%             %[obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
%             [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts);
%             obj = assign_perm(obj, pattern{1}, [0, 0, 1, 1]);
    %
%         if obj.testing == 1
%             obj = cell2matrix(obj);
%             calculate_error(obj, [
%                     0, 1, 0, 0;
%                     0, 1, 1, 0;
%                     0, 1, 1, 0;
%                     0, 0, 0, 0], struct,1)
% 
%             calculate_error(obj, [
%                     0, 1, 1, 0;
%                     0, 1, 1, 0;
%                     0, 1, 1, 0;
%                     0, 0, 0, 0], struct,1)
%             calculate_error(obj, [
%                     0, 1, 0, 0;
%                     0, 1, 1, 1;
%                     0, 1, 1, 0;
%                     0, 0, 0, 0], struct,1)
%             calculate_error(obj, [
%                     0, 0, 1, 0;
%                     1, 1, 1, 0;
%                     0, 1, 1, 1;
%                     0, 1, 0, 0], struct,1)
%         end

    %     %long extension
    %     [map, ~] = create_map([
    %                         0, 0, 1, 0;
    %                         0, 0, 1, 0;
    %                         0, 0, 1, 1;
    %                         0, 0, 1, 1]);
    %     pattern = {[0, 2, a, a]};
    %     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %     obj = assign_perm(obj, pattern{1},[0,0,1,1]);
    %
    %     [map, ~] = create_map([
    %                         0, 0, 1, 0;
    %                         0, 0, 1, 0;
    %                         0, 1, 1, 1;
    %                         0, 0, 1, 1]);
    %     pattern = {[1, 2, a, a]};
    %     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %     obj = assign_perm(obj, pattern{1},[0,0,1,1]);
    % %
    %     [map, ~] = create_map([
    %                         0, 0, 1, 0;
    %                         0, 0, 1, 0;
    %                         1, 1, 1, 1;
    %                         0, 0, 1, 1;
    %                         ]);
    %     pattern = {[2, 2, a, a]};
    %     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %     obj = assign_perm(obj, pattern{1},[0,0,1,1]);
    %
    %     if obj.testing == 1
    %         calculate_error(obj, [
    %                         0, 0, 1, 0;
    %                         0, 0, 1, 0;
    %                         0, 0, 1, 1;
    %                         0, 0, 1, 1;
    %                         ], struct)
    %     end

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


function error_code = check_lin_error(obj)
    err = calculate_error(obj, 1:6, struct("numbered",true, "h_cyclic",1   ));
    error_code =  err > 1+1e-3;
    
end


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
%
% function test_rot_180(A, B)
%     A = obj.PEPO_cell{A(1) + 1, A(2) + 1, A(3) + 1, A(4) + 1};
%     B = obj.PEPO_cell{B(1) + 1, B(2) + 1, B(3) + 1, B(4) + 1};
%     err = svds(reshape(B - permute(A, [1, 2, 5, 6, 3, 4]), [], 1))
% end
