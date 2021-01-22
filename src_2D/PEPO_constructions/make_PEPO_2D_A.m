function obj = make_PEPO_2D_A(obj)
    d = obj.dim;
    ln_prefact = obj.nf;

    %%%%%%%%%%%%%% LEVEL 1 %%%%%%%%%%%%%%
    obj.current_max_index = 1;
    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^2];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_horiz, d^2];

    %0--|--1--|--0 and all other veriants
    n=2;
    [map, ~] = create_map(1:n, obj.numopts);
    pattern = {[0,0,1,0],[1,0,0,0]};

    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    %copy to vertical cells
    obj.PEPO_cell{1, 2, 1, 1} = reshape(  obj.PEPO_cell{2, 1, 1, 1}  , [d, d, 1,  d^2,1, 1]); %right
    obj.PEPO_cell{1, 1, 1, 2} = reshape(  obj.PEPO_cell{1, 1, 2, 1}  , [d, d, 1,  1,1, d^2]); %right

    if obj.testing == 1
        err = calculate_error(obj, 1:n, obj.numopts)
        err = calculate_error(obj, (1:n)', obj.numopts)
    end

    % 0--|--1--|1--|--0 and all other veriants
    n=3;
    [map, ~] = create_map(1:n, obj.numopts);
    pattern = {[1,0,1,0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    block_11 = obj.PEPO_cell{2,1,2,1};

    %copy to equivalent blocks
    obj.PEPO_cell{2, 1, 1, 2} = reshape(block_11, [d, d, d^2, 1, 1, d^2]);
    obj.PEPO_cell{1, 2, 2, 1} = reshape(block_11, [d, d, 1, d^2, d^2, 1]);
    obj.PEPO_cell{1, 2, 1, 2} = reshape(block_11, [d, d, 1, d^2, 1, d^2]);


    if obj.testing == 1
        calculate_error(obj,[1 2 3], obj.numopts)
        calculate_error(obj,[1 2; 0 3], obj.numopts)
        calculate_error(obj,[1 0; 2 3], obj.numopts)
        calculate_error(obj,[1; 2; 3; ], obj.numopts)
    end

    %special 1
    [map, ~] = create_map([1,2;
                           3,0], obj.numopts);
    pattern = {[0,0,1,1]};                       
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    %special 2
    [map, ~] = create_map([0,2;
                           3,1], obj.numopts);
    pattern = {[1,1,0,0]}; 
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj,[1,2;3,0], obj.numopts)
        calculate_error(obj,[0,2; 3,1], obj.numopts)
    end

    %% 1 block with 3 external legs
    [map, ~] = create_map([0,2,0;
                           3,1,4], obj.numopts);
    pattern = {[1,1,1,0]}; 
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([3,2,4;
                           0,1,0], obj.numopts);
    pattern = {[1,0,1,1]}; 
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([1,0;
                           2,3;
                           4,0], obj.numopts);
    pattern = {[0,1,1,1]}; 
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0,1;
                           2,3;
                           0,4], obj.numopts);
    pattern = {[1,1,0,1]}; 
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj,[0,2,0;3,1,4], obj.numopts)
        calculate_error(obj,[3,2,4;0,1,0], obj.numopts)
        calculate_error(obj,[1,0;2,3;4,0], obj.numopts)
        calculate_error(obj,[0,1;2,3;0,4], obj.numopts)
    end
    %% 1 block with 4 external legs
    [map, ~] = create_map([0,1,0;
                           2,3,5;
                           0,4,0], obj.numopts);
    pattern = {[1,1,1,1]}; 
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    if obj.testing == 1
        calculate_error(obj,[0,1,0;2,3,5;0,4,0], obj.numopts)
    end

    %%%%%%%%%%%%%% LEVEL 2 %%%%%%%%%%%%%%
    % obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^4];
    % obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_horiz, d^4];
    % obj.current_max_index = 2;
    % %--1--|--2--|--1 and variants
    % [map, ~] = create_map(1:4, obj.numopts);
    % pattern = {[1,0,2,0],[2,0,1,0]};
    % [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    % %copy to vertical cells
    % obj.PEPO_cell{1, 3, 1, 2} = reshape(  obj.PEPO_cell{3, 1, 2, 1}  , [d, d, 1,  d^4,1, d^2]); %right
    % obj.PEPO_cell{1, 2, 1, 3} = reshape(  obj.PEPO_cell{2, 1, 3, 1}  , [d, d, 1,  d^2,1, d^4]); %right

    % if obj.testing == 1
    %     err = calculate_error(obj, 1:4, obj.numopts)
    %     err = calculate_error(obj, (1:4)', obj.numopts)
    % end
    % %--2--|--2-- and variants

    % [map, ~] = create_map(1:5, obj.numopts);
    % pattern = {[2,0,2,0]};
    % [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    % block_22 = obj.PEPO_cell{3,1,3,1};

    % %copy to equivalent blocks
    % obj.PEPO_cell{3, 1, 1, 3} = reshape(block_22, [d, d, d^4, 1, 1, d^4]);
    % obj.PEPO_cell{1, 3, 3, 1} = reshape(block_22, [d, d, 1, d^4, d^4, 1]);
    % obj.PEPO_cell{1, 3, 1, 3} = reshape(block_22, [d, d, 1, d^4, 1, d^4]);




    % if obj.testing == 1
    %     calculate_error(obj,[1 2 3;0 0 4; 0 0 5], obj.numopts) 
    %     calculate_error(obj,[1 2 3 4 5], obj.numopts)
    %     calculate_error(obj,[1 2 3 0;0 0 4 5], obj.numopts)
    %     calculate_error(obj,[2,3,4,1;0,0,0,5], obj.numopts)
    %     calculate_error(obj,[2,3,4,1;5,0,0,0], obj.numopts)
    % end

    % %special 1
    % [map, ~] = create_map([1,2;
    %                        3,0], obj.numopts);
    % pattern = {[0,0,1,1]};                       
    % [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    % %special 2
    % [map, ~] = create_map([0,2;
    %                        3,1], obj.numopts);
    % pattern = {[1,1,0,0]}; 
    % [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    % if obj.testing == 1
    %     calculate_error(obj,[1,2;3,0], obj.numopts)
    %     calculate_error(obj,[0,2; 3,1], obj.numopts)
    % end


    %%%%%%%%%%%%%% loops %%%%%%%%%%%%%%
    obj.current_max_index = obj.current_max_index+1;
    level =obj.current_max_index;
    loop_dim = d^2+1;
    %obj.cycle_index = level;

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, loop_dim];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_horiz, loop_dim];

    [map, ~] = create_map([4,2;3,1], obj.numopts);

    obj.PEPO_cell{level+1,level+1,1,1}=rand([d,d,loop_dim,loop_dim,1,1]);
    obj.PEPO_cell{level+1,1,1,level+1}=rand([d,d,loop_dim,1,1,loop_dim]);
    obj.PEPO_cell{1,level+1,level+1,1}=rand([d,d,1,loop_dim,loop_dim,1]);
    obj.PEPO_cell{1,1,level+1,level+1}=rand([d,d,1,1,loop_dim,loop_dim]);

    % 4 corner pieces
    pattern = {[level,level,0,0],[level,0,0,level],[0,level,level,0],[0,0,level,level] };

    [target, ln_prefact] = H_exp(obj, map, ln_prefact, true);
    target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(d^2, map.N));
    target = reshape(target, [d^map.N, d^map.N]);

    mul_factor = exp(ln_prefact - obj.nf);

    con_cells = get_valid_contractions(obj, map, struct('max_index', obj.current_max_index, 'pattern', {pattern}));

    x_cell = solve_non_lin(obj, pattern, {map}, {target_site}, {con_cells}, struct(), ln_prefact  )

    obj.PEPO_cell{level+1,level+1,1,1}= x_cell{1};
    obj.PEPO_cell{level+1,1,1,level+1}= x_cell{2};
    obj.PEPO_cell{1,level+1,level+1,1}= x_cell{3};
    obj.PEPO_cell{1,1,level+1,level+1}= x_cell{4};

    if obj.testing == 1
        calculate_error(obj,[1,2;3,4], obj.numopts)
    end

end
