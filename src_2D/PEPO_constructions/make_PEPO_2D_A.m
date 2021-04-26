function [obj, error_code] = make_PEPO_2D_A(obj)
    d = obj.dim;
    ln_prefact = obj.nf;

    error_code = 0;

    rot_180 = {{[3, 4, 1, 2]}};
    nl_opts = struct('Gradient', true, 'Display', 'None');

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

    
    function [obj,ln_prefact,error_code] = add_lin(obj, pattern ,ln_prefact,nl)
        if nargin < 4
           nl=0; 
        end

        [map1, ~] = create_map(  make_cross( pattern )  );

        if nl==0
            [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map1, {pattern}, ln_prefact);
        else
            [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map1}, {pattern}, ln_prefact, nl_opts, rot_180);
        end

        obj = assign_perm(obj, pattern);
        error_code = check_lin_error(obj);


        if obj.testing == 1
            err1 = calculate_error(obj, map1);
            if err1 > obj.err_tol
                disp(err1);
            end
        end  
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Blocks 1D chain %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for l = 1:2
        obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^(2*l)];
        obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^(2*l)];
        obj.current_max_index = l;
        
        
        [obj,ln_prefact,~] = add_lin(obj, [l-1,0,l,0] ,ln_prefact,1); 
        if l==1 && error_code ==1 %try again with complex numbers
            obj.complex = true;
            [obj,ln_prefact,error_code] = add_lin(obj, [l-1,0,l,0] ,ln_prefact,1);
        end
        
        [obj,ln_prefact,~] = add_lin(obj, [l,l,0,0] ,ln_prefact); 
        
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Block with 3/4 legs %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    obj.current_max_index = obj.max_index;


    %level 1
    [obj,ln_prefact,~] = add_lin(obj, [1,1,1,0] ,ln_prefact); 
    [obj,ln_prefact,~] = add_lin(obj, [1,1,1,1] ,ln_prefact); 
    
    %level 2
    [obj,ln_prefact,~] = add_lin(obj, [2,1,1,0] ,ln_prefact);
    [obj,ln_prefact,~] = add_lin(obj, [2,1,1,1] ,ln_prefact);
    
    [obj,ln_prefact,~] = add_lin(obj, [2,2,1,0] ,ln_prefact);
    [obj,ln_prefact,~] = add_lin(obj, [2,2,1,1] ,ln_prefact);
    
    [obj,ln_prefact,~] = add_lin(obj, [2,2,2,0] ,ln_prefact);
    [obj,ln_prefact,~] = add_lin(obj, [2,2,2,1] ,ln_prefact);
    [obj,ln_prefact,~] = add_lin(obj, [2,2,2,2] ,ln_prefact);


    %%
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%% LOOPS %%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %x = 2;
%     x = 2;
%     a = 3;
%     b = a + 1;
% 
%     %obj.PEPO_cell{zero_level+1,zero_level+1,zero_level+1,zero_level+1}= obj.PEPO_cell{1,1,1,1};
% 
%     alpha_dim = 8;
%     %beta_dim = 20;
% 
%     %obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^2, alpha_dim];
%     %obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^2, alpha_dim];
% 
%     obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, alpha_dim];
%     obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, alpha_dim];
% 
%     obj.current_max_index = numel(obj.virtual_level_sizes_horiz);
%     obj.max_index = obj.current_max_index;
% 
%     %simple loop
% 
%     [map, ~] = create_map([1, 2; 3, 4], obj.numopts);
% 
%     pattern = {[a, a, 0, 0]};
%     alpha_pattern_perm = {{...
%                             [2, 3, 4, 1], ...
%                             [3, 4, 1, 2], ...
%                             [4, 1, 2, 3]}};
% 
%     %nl_opts = struct('Gradient', true, 'Display', 'iter-detailed');
% 
%     lnopts = struct('Display', 0, 'maxit', 1);
%     [obj, ln_prefact] = solve_non_lin_and_assign(obj, {map}, pattern, ln_prefact, nl_opts, alpha_pattern_perm);
%     %obj = assign_perm(obj, alpha_pattern{1},a);
% 
%     if obj.testing == 1
%         err = calculate_error(obj, map);
%         if err > obj.err_tol
%             disp(err);
%         end
%     end
% 
%     obj = cell2matrix(obj);
%     err = calculate_error(obj, [
%                             1, 1, 1;
%                             1, 1, 1; ], struct, 1);
% 
%     if err > 1
%         error_code = 1;
%     end

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
    err = calculate_error(obj, 1:6, struct("numbered", true, "h_cyclic", 1));
    error_code = err > 1 + 1e-3;
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


