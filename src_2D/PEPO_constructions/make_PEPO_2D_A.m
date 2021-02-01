function obj = make_PEPO_2D_A(obj)
    d = obj.dim;
    ln_prefact = obj.nf;

    %%%%%%%%%%%%%% LEVEL 1 %%%%%%%%%%%%%%
    obj.current_max_index = 1;
    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^2];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_horiz, d^2];

    %0--|--1--|--0 and all other veriants
    n = 2;
    [map, ~] = create_map(1:n, obj.numopts);
    pattern = {[0, 0, 1, 0], [1, 0, 0, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    obj.PEPO_cell{1, 2, 1, 1} = reshape(obj.PEPO_cell{2, 1, 1, 1}, [d, d, 1, d^2, 1, 1]); %right
    obj.PEPO_cell{1, 1, 1, 2} = reshape(obj.PEPO_cell{1, 1, 2, 1}, [d, d, 1, 1, 1, d^2]); %right

    if obj.testing == 1
        err = calculate_error(obj, 1:n, obj.numopts)
        err = calculate_error(obj, (1:n)', obj.numopts)
    end

    % 0--|--1--|1--|--0 and all other veriants
    n = 3;
    [map, ~] = create_map(1:n, obj.numopts);
    pattern = {[1, 0, 1, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    block_11 = obj.PEPO_cell{2, 1, 2, 1};

    %copy to equivalent blocks
    obj.PEPO_cell{2, 1, 1, 2} = reshape(block_11, [d, d, d^2, 1, 1, d^2]);
    obj.PEPO_cell{1, 2, 2, 1} = reshape(block_11, [d, d, 1, d^2, d^2, 1]);
    obj.PEPO_cell{1, 2, 1, 2} = reshape(block_11, [d, d, 1, d^2, 1, d^2]);

    if obj.testing == 1
        calculate_error(obj, [1 2 3], obj.numopts)
        calculate_error(obj, [1 2; 0 3], obj.numopts)
        calculate_error(obj, [1 0; 2 3], obj.numopts)
        calculate_error(obj, [1; 2; 3; ], obj.numopts)
    end

    %special 1
    [map, ~] = create_map([1, 2;
                        3, 0], obj.numopts);
    pattern = {[0, 0, 1, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    %special 2
    [map, ~] = create_map([0, 2;
                        3, 1], obj.numopts);
    pattern = {[1, 1, 0, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [1, 2; 3, 0], obj.numopts)
        calculate_error(obj, [0, 2; 3, 1], obj.numopts)
    end

    %% 1 block with 3 external legs
    [map, ~] = create_map([0, 2, 0;
                        3, 1, 4], obj.numopts);
    pattern = {[1, 1, 1, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([3, 2, 4;
                        0, 1, 0], obj.numopts);
    pattern = {[1, 0, 1, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([1, 0;
                        2, 3;
                        4, 0], obj.numopts);
    pattern = {[0, 1, 1, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1;
                        2, 3;
                        0, 4], obj.numopts);
    pattern = {[1, 1, 0, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [0, 2, 0; 3, 1, 4], obj.numopts)
        calculate_error(obj, [3, 2, 4; 0, 1, 0], obj.numopts)
        calculate_error(obj, [1, 0; 2, 3; 4, 0], obj.numopts)
        calculate_error(obj, [0, 1; 2, 3; 0, 4], obj.numopts)
    end
    %% 1 block with 4 external legs
    [map, ~] = create_map([0, 1, 0;
                        2, 3, 5;
                        0, 4, 0], obj.numopts);
    pattern = {[1, 1, 1, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    if obj.testing == 1
        calculate_error(obj, [0, 1, 0; 2, 3, 5; 0, 4, 0], obj.numopts)
    end
    %%
    %%%%%%%%%%%%%% LEVEL 2 %%%%%%%%%%%%%%
    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^4];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_horiz, d^4];
    obj.current_max_index = 2;
    %%
    %%%%--1--|--2--|--1 and variants%%%%
    [map, ~] = create_map(1:4, obj.numopts);
    pattern = {[1, 0, 2, 0], [2, 0, 1, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    %%horizontal%%
    %equivalent
    obj.PEPO_cell{3, 1, 1, 2} = reshape(obj.PEPO_cell{3, 1, 2, 1}, [d, d, d^4, 1, 1, d^2]);
    obj.PEPO_cell{1, 2, 3, 1} = reshape(obj.PEPO_cell{2, 1, 3, 1}, [d, d, 1, d^2, d^4, 1]);
    %inequivalent
    [map, ~] = create_map([1, 2, 4;
                        3, 0, 0], obj.numopts);
    pattern = {[0, 0, 2, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 4;
                        1, 2, 3], obj.numopts);
    pattern = {[2, 1, 0, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        err = calculate_error(obj, 1:4, obj.numopts)
        err = calculate_error(obj, [0, 0, 4; 1, 2, 3], obj.numopts)
        err = calculate_error(obj, [1, 2, 4; 3, 0, 0], obj.numopts)
    end

    %%vertical%%
    %copy regular sites horizontal
    obj.PEPO_cell{1, 3, 1, 2} = reshape(obj.PEPO_cell{3, 1, 2, 1}, [d, d, 1, d^4, 1, d^2]);
    obj.PEPO_cell{1, 2, 1, 3} = reshape(obj.PEPO_cell{2, 1, 3, 1}, [d, d, 1, d^2, 1, d^4]);
    %equivalent vertical
    obj.PEPO_cell{1, 3, 2, 1} = reshape(obj.PEPO_cell{1, 3, 1, 2}, [d, d, 1, d^4, d^2, 1]);
    obj.PEPO_cell{2, 1, 1, 3} = reshape(obj.PEPO_cell{1, 2, 1, 3}, [d, d, d^2, 1, 1, d^4]);
    %inquivalent vert
    [map, ~] = create_map([1, 2;
                        3, 0;
                        4, 0], obj.numopts);
    pattern = {[0, 0, 1, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 2;
                        0, 1;
                        3, 4], obj.numopts);
    pattern = {[1, 2, 0, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        err = calculate_error(obj, [1, 2; 3, 0; 4, 0], obj.numopts)
        err = calculate_error(obj, [0, 2; 0, 1; 3, 4], obj.numopts)
        err = calculate_error(obj, (1:4)', obj.numopts)
    end
    %%
    %%%%--2--|--2-- and variants%%%%

    [map, ~] = create_map(1:5, obj.numopts);
    pattern = {[2, 0, 2, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    block_22 = obj.PEPO_cell{3, 1, 3, 1};

    %copy to equivalent blocks
    obj.PEPO_cell{3, 1, 1, 3} = reshape(block_22, [d, d, d^4, 1, 1, d^4]);
    obj.PEPO_cell{1, 3, 3, 1} = reshape(block_22, [d, d, 1, d^4, d^4, 1]);
    obj.PEPO_cell{1, 3, 1, 3} = reshape(block_22, [d, d, 1, d^4, 1, d^4]);
    %inequivalent
    [map, ~] = create_map([0, 3, 4, 5; 1, 2, 0, 0], obj.numopts);
    pattern = {[0, 0, 2, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %inequivalent
    [map, ~] = create_map([0, 0, 4, 5; 1, 2, 3, 0], obj.numopts);
    pattern = {[2, 2, 0, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [1 2 3; 0 0 4; 0 0 5], obj.numopts)
        calculate_error(obj, [1 2 3 4 5], obj.numopts)
        calculate_error(obj, [1 2 3 0; 0 0 4 5], obj.numopts)
        calculate_error(obj, [2, 3, 4, 1; 0, 0, 0, 5], obj.numopts)
        calculate_error(obj, [2, 3, 4, 1; 5, 0, 0, 0], obj.numopts)

        calculate_error(obj, [0, 3, 4, 5; 1, 2, 0, 0], obj.numopts)
        calculate_error(obj, [0, 0, 4, 5; 1, 2, 3, 0], obj.numopts)
    end
    %%
    %1 2 and 2 1 levels
    %horiz left
    [map, ~] = create_map([2, 3, 4, 5;
                        0, 1, 0, 0], obj.numopts);
    pattern = {[1, 0, 2, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 3, 0, 0;
                        2, 1, 4, 5], obj.numopts);
    pattern = {[1, 1, 2, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([3, 0, 0;
                        2, 4, 5;
                        1, 0, 0], obj.numopts);
    pattern = {[0, 1, 2, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %horiz right
    [map, ~] = create_map([2, 3, 4, 5;
                        0, 0, 1, 0], obj.numopts);
    pattern = {[2, 0, 1, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 3, 0;
                        2, 1, 4, 5], obj.numopts);
    pattern = {[2, 1, 1, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 3;
                        2, 4, 5;
                        0, 0, 1], obj.numopts);
    pattern = {[2, 1, 0, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    %vert up
    [map, ~] = create_map([1, 0;
                        2, 3;
                        4, 0;
                        5, 0], obj.numopts);
    pattern = {[0, 1, 1, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1;
                        2, 3;
                        0, 4;
                        0, 5], obj.numopts);
    pattern = {[1, 1, 0, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([2, 3, 1;
                        0, 4, 0;
                        0, 5, 0], obj.numopts);
    pattern = {[1, 0, 1, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %vert down
    [map, ~] = create_map([1, 0;
                        2, 0;
                        4, 3;
                        5, 0], obj.numopts);
    pattern = {[0, 2, 1, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1;
                        0, 3;
                        2, 4;
                        0, 5], obj.numopts);
    pattern = {[1, 2, 0, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 3, 0;
                        0, 4, 0;
                        2, 5, 1], obj.numopts);
    pattern = {[1, 2, 1, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [0, 3, 0; 0, 4, 0; 2, 5, 1], obj.numopts)
        calculate_error(obj, [1, 0; 2, 0; 4, 3; 5, 0], obj.numopts)
    end

    %% 1 2 level and 3 1 levels
    %horiz
    [map, ~] = create_map([0, 5, 0, 0;
                        1, 2, 3, 4;
                        0, 6, 0, 0], obj.numopts);
    pattern = {[1, 1, 2, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 5, 0;
                        1, 2, 3, 4;
                        0, 0, 6, 0], obj.numopts);
    pattern = {[2, 1, 1, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %vert
    [map, ~] = create_map([0, 4, 0;
                        1, 2, 3;
                        0, 5, 0;
                        0, 6, 0], obj.numopts);
    pattern = {[1, 1, 1, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 4, 0;
                        0, 2, 0;
                        1, 5, 3;
                        0, 6, 0], obj.numopts);
    pattern = {[1, 2, 1, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [0, 4, 0; 1, 2, 3; 0, 5, 0; 0, 6, 0], obj.numopts)
    end

    %% 2 2 levels + 1 level
    [map, ~] = create_map([1, 2, 3, 4, 5;
                        0, 0, 6, 0, 0], obj.numopts);
    pattern = {[2, 0, 2, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 3, 0, 0;
                        1, 2, 6, 4, 5], obj.numopts);
    pattern = {[2, 1, 2, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([1, 0;
                        2, 0;
                        3, 4;
                        5, 0;
                        6, 0], obj.numopts);
    pattern = {[0, 2, 1, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1;
                        0, 2;
                        3, 4;
                        0, 5;
                        0, 6], obj.numopts);
    pattern = {[1, 2, 0, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [1, 2, 3, 4, 5; 0, 0, 6, 0, 0], obj.numopts)
        calculate_error(obj, [0, 1; 0, 2; 3, 4; 0, 5; 0, 6], obj.numopts)

        %         calculate_error(obj, [0, 0, 3, 0, 0, 0;
        %                         1, 2, 5, 6, 7, 8;
        %                         0, 0, 0, 4, 0, 0], obj.numopts)

    end

    %% 2 2 levels + 2 1 level
    [map, ~] = create_map([0, 0, 6, 0, 0;
                        1, 2, 3, 4, 5;
                        0, 0, 7, 0, 0], obj.numopts);
    pattern = {[2, 1, 2, 1]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1, 0;
                        0, 2, 0;
                        6, 3, 7;
                        0, 4, 0;
                        0, 5, 0], obj.numopts);
    pattern = {[1, 2, 1, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [0, 0, 6, 0, 0; 1, 2, 3, 4, 5; 0, 0, 7, 0, 0], obj.numopts)
    end

    %%
    % 3 2 levels
    [map, ~] = create_map([1, 2, 3, 4, 5;
                        0, 0, 6, 0, 0
                        0, 0, 7, 0, 0], obj.numopts);
    pattern = {[2, 0, 2, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 7, 0, 0;
                        0, 0, 3, 0, 0;
                        1, 2, 6, 4, 5], obj.numopts);
    pattern = {[2, 2, 2, 0]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([1, 0, 0;
                        2, 0, 0;
                        3, 4, 7;
                        5, 0, 0;
                        6, 0, 0], obj.numopts);
    pattern = {[0, 2, 2, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 1;
                        0, 0, 2;
                        7, 3, 4;
                        0, 0, 5;
                        0, 0, 6], obj.numopts);
    pattern = {[2, 2, 0, 2]};
    [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [1, 2, 3, 4, 5; 0, 0, 6, 0, 0; 0, 0, 7, 0, 0], obj.numopts)
        calculate_error(obj, [0, 0, 7, 0, 0; 0, 0, 3, 0, 0; 1, 2, 6, 4, 5], obj.numopts)

    end
    %% has almost no influence bc rank deficient
    %3 2 legs and 1 leg
    % [map, ~] = create_map([0, 0, 8, 0, 0;
    %                     1, 2, 3, 4, 5;
    %                     0, 0, 6, 0, 0
    %                     0, 0, 7, 0, 0], obj.numopts);
    % pattern = {[2, 1, 2, 2]};
    % [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    % %
    % [map, ~] = create_map([0, 0, 8, 0, 0;
    %                     0, 0, 3, 0, 0;
    %                     1, 2, 6, 4, 5
    %                     0, 0, 7, 0, 0], obj.numopts);
    % pattern = {[2, 2, 2, 1]};
    % [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    % %
    % [map, ~] = create_map([0, 1, 0, 0;
    %                     0, 2, 0, 0;
    %                     8, 3, 4, 7;
    %                     0, 5, 0, 0;
    %                     0, 6, 0, 0], obj.numopts);
    % pattern = {[1, 2, 2, 2]};
    % [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    % %
    % [map, ~] = create_map([0, 0, 1, 0;
    %                     0, 0, 2, 0;
    %                     7, 3, 4, 8;
    %                     0, 0, 5, 0;
    %                     0, 0, 6, 0], obj.numopts);
    % pattern = {[2, 2, 1, 2]};
    % [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    % if obj.testing == 1
    %     calculate_error(obj, [0, 0, 8, 0, 0; 1, 2, 3, 4, 5; 0, 0, 6, 0, 0; 0, 0, 7, 0, 0], obj.numopts)
    % end
    %%
    %4 leg 2 tensors
    %     [map, ~] = create_map([0,0,1,0,0;
    %                            0,0,2,0,0;
    %                            7,3,4,8,9;
    %                            0,0,5,0,0;
    %                            0,0,6,0,0], obj.numopts);
    %     pattern = {[2,2,1,2]};
    %     [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    %%
    %%%%%%%%%%%%%% loops %%%%%%%%%%%%%%
    obj.current_max_index = obj.current_max_index + 1;
    level = obj.current_max_index;
    obj.max_index = 3;
    loop_dim = d^2 + 2;
    
    int_level = 2;
    int_dim = obj.virtual_level_sizes_horiz(int_level+1);
    %obj.cycle_index = level;

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, loop_dim];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_horiz, loop_dim];

    %simple loop
    [map1, ~] = create_map([1, 2; 3, 4], obj.numopts);
    pattern1 = {[level, level, 0, 0], [level, 0, 0, level], [0, level, level, 0], [0, 0, level, level]};

    [obj, ln_prefact] = solve_non_lin_and_assign(obj, map1, pattern1, ln_prefact, loop_dim, 1);

    if obj.testing == 1
        calculate_error(obj, [1, 2; 3, 4], obj.numopts)
    end

%     %6 loop horizontal
%     [map2, ~] = create_map([1, 4, 5; 3, 2, 6], obj.numopts);
%     pattern2 = {[level, int_level, level, 0], [level, 0, level, int_level]};
% 
%     [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map2, pattern2, ln_prefact, int_dim, 1);
% 
%     if obj.testing == 1
%         calculate_error(obj, [1, 2, 5; 3, 4, 6], obj.numopts)
%     end
%     
%    
%     
% 
%     %6 loop vertical
%     [map3, ~] = create_map([1, 2; 3, 4; 5, 6], obj.numopts);
%     pattern3 = {[0, level, int_level, level], [int_level, level, 0, level]};
% 
%     [obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map3, pattern3, ln_prefact, int_dim,1);
% 
%     if obj.testing == 1
%         calculate_error(obj, [1, 2; 3, 4; 5, 6], obj.numopts)
%     end

    % solve simulataneous linear problem with simple loop and 2 double
    % loops
    
    

%     [target1, ~] = H_exp(obj, map1, ln_prefact, false);
%     target_site1 = reshape(permute(target1, site_ordering_permute(map1.N)), dimension_vector(d^2, map1.N));
%     
%     [target2, ~] = H_exp(obj, map2, ln_prefact, false);
%     target_site2 = reshape(permute(target2, site_ordering_permute(map2.N)), dimension_vector(d^2, map2.N));
%     
%     [target3, ~] = H_exp(obj, map2, ln_prefact, false);
%     target_site3 = reshape(permute(target3, site_ordering_permute(map3.N)), dimension_vector(d^2, map3.N));
%     
%     con_cells1 = get_valid_contractions(obj, map1, struct('max_index', obj.current_max_index));
%     con_cells2 = get_valid_contractions(obj, map2, struct('max_index', obj.current_max_index));
%     con_cells3 = get_valid_contractions(obj, map3, struct('max_index', obj.current_max_index));
%     
%     patterns = [pattern1,pattern2,pattern3];
%     maps = {map1, map2, map3};
%     con_cells = {con_cells1,con_cells2,con_cells3};
%     targets = {target_site1,target_site2,target_site3};
%     
%      x_cell = solve_non_lin(obj, patterns, maps, targets, con_cells, struct(), ln_prefact)
    
    %2by2 square
    %      [map, ~] = create_map([1,2,3;4,5,6;7,8,9], obj.numopts);
    %      pattern = {[level,level,level,level] };
    %
    %      %[obj, ~, ~, ln_prefact, rank_x] = solve_lin_and_assign(obj, map, pattern, ln_prefact,loop_dim);
    %      [obj, ln_prefact]= solve_non_lin_and_assign(obj, map, pattern, ln_prefact,loop_dim);
    %
    %

    %      if obj.testing == 1
    %          calculate_error(obj,[1,2,3;4,5,6;7,8,9], obj.numopts)
    %      end

end
