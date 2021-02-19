function obj = make_PEPO_2D_A(obj)
    d = obj.dim;
    ln_prefact = obj.nf;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Block with 1/2 legs %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%% LEVEL 1 %%%%%%%%%%%%%%
    obj.current_max_index = 1;
    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^2];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^2];

    %0--|--1--|--0 and all other veriants
    n = 2;
    [map, ~] = create_map(1:n, obj.numopts);
    pattern = {[0, 0, 1, 0], [1, 0, 0, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    obj.PEPO_cell{1, 2, 1, 1} = reshape(obj.PEPO_cell{2, 1, 1, 1}, [d, d, 1, d^2, 1, 1]); %right
    obj.PEPO_cell{1, 1, 1, 2} = reshape(obj.PEPO_cell{1, 1, 2, 1}, [d, d, 1, 1, 1, d^2]); %right

    if obj.testing == 1
        calculate_error(obj, 1:n, obj.numopts)
        calculate_error(obj, (1:n)', obj.numopts)
    end

    % 0--|--1--|1--|--0 and all other veriants
    %#6
    n = 3;
    [map, ~] = create_map(1:n, obj.numopts);
    pattern = {[1, 0, 1, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

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
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    %special 2
    [map, ~] = create_map([0, 2;
                        3, 1], obj.numopts);
    pattern = {[1, 1, 0, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [1, 2; 3, 0], obj.numopts)
        calculate_error(obj, [0, 2; 3, 1], obj.numopts)
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
    pattern = {[1, 0, 2, 0], [2, 0, 1, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    %%horizontal%%
    %equivalent
    obj.PEPO_cell{3, 1, 1, 2} = reshape(obj.PEPO_cell{3, 1, 2, 1}, [d, d, d^4, 1, 1, d^2]);
    obj.PEPO_cell{1, 2, 3, 1} = reshape(obj.PEPO_cell{2, 1, 3, 1}, [d, d, 1, d^2, d^4, 1]);
    %inequivalent
    [map, ~] = create_map([1, 2, 4;
                        3, 0, 0], obj.numopts);
    pattern = {[0, 0, 2, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 4;
                        1, 2, 3], obj.numopts);
    pattern = {[2, 1, 0, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, 1:4, obj.numopts)
        calculate_error(obj, [0, 0, 4; 1, 2, 3], obj.numopts)
        calculate_error(obj, [1, 2, 4; 3, 0, 0], obj.numopts)
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
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 2;
                        0, 1;
                        3, 4], obj.numopts);
    pattern = {[1, 2, 0, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [1, 2; 3, 0; 4, 0], obj.numopts)
        calculate_error(obj, [0, 2; 0, 1; 3, 4], obj.numopts)
        calculate_error(obj, (1:4)', obj.numopts)
    end
    %%
    %%%%--2--|--2-- and variants%%%%
    %# 4*3/2= 6

    [map, ~] = create_map(1:5, obj.numopts);
    pattern = {[2, 0, 2, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    block_22 = obj.PEPO_cell{3, 1, 3, 1};

    %copy to equivalent blocks
    obj.PEPO_cell{3, 1, 1, 3} = reshape(block_22, [d, d, d^4, 1, 1, d^4]);
    obj.PEPO_cell{1, 3, 3, 1} = reshape(block_22, [d, d, 1, d^4, d^4, 1]);
    obj.PEPO_cell{1, 3, 1, 3} = reshape(block_22, [d, d, 1, d^4, 1, d^4]);
    %inequivalent
    [map, ~] = create_map([0, 3, 4, 5; 1, 2, 0, 0], obj.numopts);
    pattern = {[0, 0, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %inequivalent
    [map, ~] = create_map([0, 0, 4, 5; 1, 2, 3, 0], obj.numopts);
    pattern = {[2, 2, 0, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [1 2 3; 0 0 4; 0 0 5], obj.numopts)
        calculate_error(obj, [1 2 3 4 5], obj.numopts)
        calculate_error(obj, [1 2 3 0; 0 0 4 5], obj.numopts)
        calculate_error(obj, [2, 3, 4, 1; 0, 0, 0, 5], obj.numopts)
        calculate_error(obj, [2, 3, 4, 1; 5, 0, 0, 0], obj.numopts)

        calculate_error(obj, [0, 3, 4, 5; 1, 2, 0, 0], obj.numopts)
        calculate_error(obj, [0, 0, 4, 5; 1, 2, 3, 0], obj.numopts)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Block with 3/4 legs %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    obj.current_max_index = obj.max_index;

    %%%%%%%%%%%%%% LEVEL 1 %%%%%%%%%%%%%%
    %% 1 block with 3 external 1 legs
    %# 4
    [map, ~] = create_map([0, 2, 0;
                        3, 1, 4], obj.numopts);
    pattern = {[1, 1, 1, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([3, 2, 4;
                        0, 1, 0], obj.numopts);
    pattern = {[1, 0, 1, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([1, 0;
                        2, 3;
                        4, 0], obj.numopts);
    pattern = {[0, 1, 1, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1;
                        2, 3;
                        0, 4], obj.numopts);
    pattern = {[1, 1, 0, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

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
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    if obj.testing == 1
        calculate_error(obj, [0, 1, 0; 2, 3, 5; 0, 4, 0], obj.numopts)
    end

    %%%%%%%%%%%%%% LEVEL 2 %%%%%%%%%%%%%%

    %%
    %1 2 and 2 1 levels
    %#12
    %horiz left
    [map, ~] = create_map([2, 3, 4, 5;
                        0, 1, 0, 0], obj.numopts);
    pattern = {[1, 0, 2, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 3, 0, 0;
                        2, 1, 4, 5], obj.numopts);
    pattern = {[1, 1, 2, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([3, 0, 0;
                        2, 4, 5;
                        1, 0, 0], obj.numopts);
    pattern = {[0, 1, 2, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %horiz right
    [map, ~] = create_map([2, 3, 4, 5;
                        0, 0, 1, 0], obj.numopts);
    pattern = {[2, 0, 1, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 3, 0;
                        2, 1, 4, 5], obj.numopts);
    pattern = {[2, 1, 1, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 3;
                        2, 4, 5;
                        0, 0, 1], obj.numopts);
    pattern = {[2, 1, 0, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    %vert up
    [map, ~] = create_map([1, 0;
                        2, 3;
                        4, 0;
                        5, 0], obj.numopts);
    pattern = {[0, 1, 1, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1;
                        2, 3;
                        0, 4;
                        0, 5], obj.numopts);
    pattern = {[1, 1, 0, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([2, 3, 1;
                        0, 4, 0;
                        0, 5, 0], obj.numopts);
    pattern = {[1, 0, 1, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %vert down
    [map, ~] = create_map([1, 0;
                        2, 0;
                        4, 3;
                        5, 0], obj.numopts);
    pattern = {[0, 2, 1, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1;
                        0, 3;
                        2, 4;
                        0, 5], obj.numopts);
    pattern = {[1, 2, 0, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 3, 0;
                        0, 4, 0;
                        2, 5, 1], obj.numopts);
    pattern = {[1, 2, 1, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [0, 3, 0; 0, 4, 0; 2, 5, 1], obj.numopts)
        calculate_error(obj, [1, 0; 2, 0; 4, 3; 5, 0], obj.numopts)
    end

    %% 1 2 loop_level and 3 1 levels
    %4
    %horiz
    [map, ~] = create_map([0, 5, 0, 0;
                        1, 2, 3, 4;
                        0, 6, 0, 0], obj.numopts);
    pattern = {[1, 1, 2, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 5, 0;
                        1, 2, 3, 4;
                        0, 0, 6, 0], obj.numopts);
    pattern = {[2, 1, 1, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %vert
    [map, ~] = create_map([0, 4, 0;
                        1, 2, 3;
                        0, 5, 0;
                        0, 6, 0], obj.numopts);
    pattern = {[1, 1, 1, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 4, 0;
                        0, 2, 0;
                        1, 5, 3;
                        0, 6, 0], obj.numopts);
    pattern = {[1, 2, 1, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [0, 4, 0; 1, 2, 3; 0, 5, 0; 0, 6, 0], obj.numopts)
    end

    %% 2 2 levels + 1 loop_level
    %#12
    [map, ~] = create_map([1, 2, 3, 4, 5;
                        0, 0, 6, 0, 0], obj.numopts);
    pattern = {[2, 0, 2, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 3, 0, 0;
                        1, 2, 6, 4, 5], obj.numopts);
    pattern = {[2, 1, 2, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([1, 0;
                        2, 0;
                        3, 4;
                        5, 0;
                        6, 0], obj.numopts);
    pattern = {[0, 2, 1, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1;
                        0, 2;
                        3, 4;
                        0, 5;
                        0, 6], obj.numopts);
    pattern = {[1, 2, 0, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    % onder rechts
    [map, ~] = create_map([4, 0, 0;
                        3, 1, 2;
                        5, 0, 0;
                        6, 0, 0], obj.numopts);
    pattern = {[0, 1, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([4, 3, 1, 2;
                        0, 5, 0, 0;
                        0, 6, 0, 0], obj.numopts);
    pattern = {[1, 0, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %onder links
    [map, ~] = create_map([0, 0, 4;
                        3, 1, 2;
                        0, 0, 5;
                        0, 0, 6], obj.numopts);
    pattern = {[2, 1, 0, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([4, 3, 1, 2;
                        0, 0, 5, 0;
                        0, 0, 6, 0], obj.numopts);
    pattern = {[2, 0, 1, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %boven links
    [map, ~] = create_map([0, 0, 4;
                        0, 0, 2;
                        3, 1, 5;
                        0, 0, 6], obj.numopts);
    pattern = {[2, 2, 0, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 4, 0;
                        0, 0, 2, 0;
                        3, 1, 5, 6; ], obj.numopts);
    pattern = {[2, 2, 1, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %boven rechts
    [map, ~] = create_map([0, 6, 0, 0;
                        0, 5, 0, 0;
                        4, 3, 1, 2; ], obj.numopts);
    pattern = {[1, 2, 2, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([6, 0, 0;
                        5, 0, 0;
                        3, 1, 2;
                        4, 0, 0; ], obj.numopts);
    pattern = {[0, 2, 2, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [1, 2, 3, 4, 5; 0, 0, 6, 0, 0], obj.numopts)
        calculate_error(obj, [0, 1; 0, 2; 3, 4; 0, 5; 0, 6], obj.numopts)

        %         calculate_error(obj, [0, 0, 3, 0, 0, 0;
        %                         1, 2, 5, 6, 7, 8;
        %                         0, 0, 0, 4, 0, 0], obj.numopts)

    end

    %% 2 2 levels + 2 1 loop_level
    %3 4*3/2=6
    [map, ~] = create_map([0, 0, 6, 0, 0;
                        1, 2, 3, 4, 5;
                        0, 0, 7, 0, 0], obj.numopts);
    pattern = {[2, 1, 2, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 5, 0;
                        0, 0, 6, 0;
                        1, 2, 3, 4;
                        0, 0, 7, 0], obj.numopts);
    pattern = {[2, 2, 1, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 5, 0, 0;
                        0, 6, 0, 0;
                        1, 2, 3, 4;
                        0, 7, 0, 0], obj.numopts);
    pattern = {[1, 2, 2, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 5, 0;
                        1, 2, 6, 4;
                        0, 0, 3, 0;
                        0, 0, 7, 0], obj.numopts);
    pattern = {[2, 1, 1, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 5, 0, 0;
                        1, 2, 6, 4;
                        0, 3, 0, 0;
                        0, 7, 0, 0], obj.numopts);
    pattern = {[1, 1, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1, 0;
                        0, 2, 0;
                        6, 3, 7;
                        0, 4, 0;
                        0, 5, 0], obj.numopts);
    pattern = {[1, 2, 1, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [0, 0, 6, 0, 0; 1, 2, 3, 4, 5; 0, 0, 7, 0, 0], obj.numopts)
    end

    %%
    % 3 2 levels
    [map, ~] = create_map([1, 2, 3, 4, 5;
                        0, 0, 6, 0, 0
                        0, 0, 7, 0, 0], obj.numopts);
    pattern = {[2, 0, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 7, 0, 0;
                        0, 0, 3, 0, 0;
                        1, 2, 6, 4, 5], obj.numopts);
    pattern = {[2, 2, 2, 0]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([1, 0, 0;
                        2, 0, 0;
                        3, 4, 7;
                        5, 0, 0;
                        6, 0, 0], obj.numopts);
    pattern = {[0, 2, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 1;
                        0, 0, 2;
                        7, 3, 4;
                        0, 0, 5;
                        0, 0, 6], obj.numopts);
    pattern = {[2, 2, 0, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [1, 2, 3, 4, 5; 0, 0, 6, 0, 0; 0, 0, 7, 0, 0], obj.numopts)
        calculate_error(obj, [0, 0, 7, 0, 0; 0, 0, 3, 0, 0; 1, 2, 6, 4, 5], obj.numopts)

    end

%     %% 3 2 legs and 1 leg %works but takes a few second
%     [map, ~] = create_map([0, 0, 8, 0, 0;
%                         1, 2, 3, 4, 5;
%                         0, 0, 6, 0, 0
%                         0, 0, 7, 0, 0], obj.numopts);
%     pattern = {[2, 1, 2, 2]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
%     %
%     [map, ~] = create_map([0, 0, 8, 0, 0;
%                         0, 0, 3, 0, 0;
%                         1, 2, 6, 4, 5
%                         0, 0, 7, 0, 0], obj.numopts);
%     pattern = {[2, 2, 2, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
%     %
%     [map, ~] = create_map([0, 1, 0, 0;
%                         0, 2, 0, 0;
%                         8, 3, 4, 7;
%                         0, 5, 0, 0;
%                         0, 6, 0, 0], obj.numopts);
%     pattern = {[1, 2, 2, 2]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
%     %
%     [map, ~] = create_map([0, 0, 1, 0;
%                         0, 0, 2, 0;
%                         7, 3, 4, 8;
%                         0, 0, 5, 0;
%                         0, 0, 6, 0], obj.numopts);
%     pattern = {[2, 2, 1, 2]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
% 
%     if obj.testing == 1
%         calculate_error(obj, [0, 0, 8, 0, 0; 1, 2, 3, 4, 5; 0, 0, 6, 0, 0; 0, 0, 7, 0, 0], obj.numopts)
%     end
    %%
    % 4 2 leg requires to much ram
    %     [map, ~] = create_map([0, 0, 1, 0, 0;
    %                         0, 0, 2, 0, 0;
    %                         7, 3, 4, 8, 9;
    %                         0, 0, 5, 0, 0;
    %                         0, 0, 6, 0, 0], obj.numopts);
    %     pattern = {[2, 2, 2 , 2 ]};
    %     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    %     if obj.testing == 1
    %         calculate_error(obj, [0, 0, 8, 0, 0; 1, 2, 3, 4, 5; 0, 0, 6, 0, 0; 0, 0, 7, 0, 0], obj.numopts)
    %     end
    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% LOOPS %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    loop_level = obj.current_max_index + 1;
    beta_level = loop_level + 1;
    gamma_level = beta_level + 2;
    obj.current_max_index = obj.current_max_index + 5;
    obj.max_index = obj.current_max_index;
    loop_dim = 11;
    beta_dim = 8;
    gamma_dim = 16;

    lnlopts = struct('Display',1,'maxit', 5);
    %obj.cycle_index = loop_level;

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, loop_dim, beta_dim, beta_dim, gamma_dim,gamma_dim];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, loop_dim, beta_dim, beta_dim, gamma_dim,gamma_dim];

    %simple loop
    [map1, ~] = create_map([1, 2; 3, 4], obj.numopts);
    pattern1 = {[loop_level, loop_level, 0, 0], [loop_level, 0, 0, loop_level], [0, loop_level, loop_level, 0], [0, 0, loop_level, loop_level]};

    %[obj, ln_prefact] = solve_non_lin_and_assign(obj, map1, pattern1, ln_prefact, struct('Algoritm', 'Levenberg-Marquardt'));
    obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, lnlopts);
    
    
    if obj.testing == 1
        calculate_error(obj, [1, 2; 3, 4], obj.numopts)
    end


    
    %simple loop
    [map1, ~] = create_map([1, 2,5; 3, 4,6], obj.numopts);
    pattern1 = {[loop_level, loop_level, loop_level, 0], [loop_level, 0, loop_level, loop_level]};

    %[obj, ln_prefact] = solve_non_lin_and_assign(obj, map1, pattern1, ln_prefact, struct('Algoritm', 'Levenberg-Marquardt'));
    obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, lnlopts);
    
    
    if obj.testing == 1
        calculate_error(obj, [1, 2,5; 3, 4,6], obj.numopts)
    end
    
    
    
    
    
    %
    fprintf("| ")

%     %3 legs loop
% 
%     nopts = struct('Display', 'iter-detailed', 'maxit', 10);
%     
%     %essentially non lin but lin is good enough
% 
%     %%%%%%%%%%%
%     %left upper
%     [map, ~] = create_map([1, 0;
%                         2, 3;
%                         4, 5; ], obj.numopts);
% 
%     pattern = {[0, 1, beta_level, beta_level], [0, beta_level, loop_level, 0], [beta_level, 0, 0, loop_level]};
% 
%     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
%     %[obj, ln_prefact] = solve_non_lin_and_assign(obj, map, pattern, ln_prefact, nopts);
%     %
%     pattern = {[0, 1, beta_level, beta_level]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([1, 2, 3;
%                         0, 4, 5; ], obj.numopts);
% 
%     pattern = {[1, 0, beta_level, beta_level]};
%     %[obj, ln_prefact] = solve_non_lin_and_assign(obj, map, pattern, ln_prefact, nopts);
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 6, 0;
%                         1, 2, 3;
%                         0, 4, 5; ], obj.numopts);
%     pattern = {[1, 1, beta_level, beta_level]};
%     %[obj, ln_prefact] = solve_non_lin_and_assign(obj, map, pattern, ln_prefact, nopts);
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %%%%%%%%%%%
%     %left lower
%     [map, ~] = create_map([0, 2, 3;
%                         1, 4, 5; ], obj.numopts);
%     pattern = {[1, beta_level + 1, beta_level + 1, 0], [0, 0, loop_level, beta_level + 1], [beta_level + 1, loop_level, 0, 0]};
%     
%     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact,lnlopts);
%     %[obj, ln_prefact] = solve_non_lin_and_assign(obj, map, pattern, ln_prefact, nopts);
%     %
%     pattern = {[1, beta_level + 1, beta_level + 1, 0]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([2, 3;
%                         4, 5;
%                         1, 0; ], obj.numopts);
%     pattern = {[0, beta_level + 1, beta_level + 1, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 2, 3;
%                         1, 4, 5;
%                         0, 6, 0; ], obj.numopts);
%     pattern = {[1, beta_level + 1, beta_level + 1, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %%%%%%%%%%%
%     %right upper
%     [map, ~] = create_map([2, 3, 1;
%                         4, 5, 0; ], obj.numopts);
%     pattern = {[beta_level + 1, 0, 1, beta_level + 1], [0, 0, beta_level + 1, loop_level], [loop_level, beta_level + 1, 0, 0]};
%     
%     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
%     %[obj, ln_prefact] = solve_non_lin_and_assign(obj, map, pattern, ln_prefact, nopts);
%     %
%     pattern = {[beta_level + 1, 0, 1, beta_level + 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 1;
%                         2, 3;
%                         4, 5; ], obj.numopts);
% 
%     pattern = {[beta_level + 1, 1, 0, beta_level + 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 1, 0;
%                         2, 3, 6;
%                         4, 5, 0; ], obj.numopts);
%     pattern = {[beta_level + 1, 1, 1, beta_level + 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %%%%%%%%%%%
%     %right lower
%     [map, ~] = create_map([2, 3, 0;
%                         4, 5, 1; ], obj.numopts);
%     pattern = {[beta_level, beta_level, 1, 0], [0, loop_level, beta_level, 0], [loop_level, 0, 0, beta_level]};
%     
%     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
%     %[obj, ln_prefact] = solve_non_lin_and_assign(obj, map, pattern, ln_prefact, nopts);
%     %
%     pattern = {[beta_level, beta_level, 1, 0]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([2, 3;
%                         4, 5;
%                         0, 1], obj.numopts);
%     pattern = {[beta_level, beta_level, 0, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([2, 3, 0;
%                         4, 5, 1;
%                         0, 6, 0], obj.numopts);
%     pattern = {[beta_level, beta_level, 1, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
% 
%     if obj.testing == 1
%         calculate_error(obj, [1, 2, 3; 0, 4, 5; ], obj.numopts)
%         calculate_error(obj, [1, 0; 2, 3; 4, 5; ], obj.numopts)
%         calculate_error(obj, [0, 6, 0; 1, 2, 3; 0, 4, 5; ], obj.numopts)
% 
%         calculate_error(obj, [0, 2, 3; 1, 4, 5; ], obj.numopts)
%         calculate_error(obj, [2, 3; 4, 5; 1, 0; ], obj.numopts)
%         calculate_error(obj, [0, 2, 3; 1, 4, 5; 0, 6, 0; ], obj.numopts)
% 
%         calculate_error(obj, [0, 2, 0; 0, 3, 4; 1, 5, 6], obj.numopts)
%         calculate_error(obj, [0, 2, 0; 7, 3, 4; 0, 5, 6; 0, 1, 0], obj.numopts)
% 
%         %fan
%         calculate_error(obj, [0, 0, 1, 0; 2, 3, 4, 0; 0, 5, 6, 7; 0, 8, 0, 0], obj.numopts)
% 
%         calculate_error(obj, [1, 0, ; 3, 4, ; 0, 2; ], obj.numopts)
% 
%         calculate_error(obj, [0, 2, 0, 0; 7, 3, 4, 0; 0, 5, 6, 8; 0, 1, 9, 0], obj.numopts)
% 
%         calculate_error(obj, [0, 0, 2; 1, 3, 4; 0, 5, 6], obj.numopts)
% 
%         calculate_error(obj, [0, 2, 3; 1, 4, 5; ], obj.numopts)
%         calculate_error(obj, [2, 3; 4, 5; 1, 0; ], obj.numopts)
%         %double corner
% 
%     end
%     %%
%     %do two corners at same time
%     %%%%%%%%%%
%     %%%%%%%%up
%     nopts = struct('Display', 'iter-detailed', 'maxit', 10);
%      lnlopts = struct('Display',1,'maxit', 100);
%     
%     [map, ~] = create_map([1, 2, 3, 6;
%                         0, 4, 5, 0; ], obj.numopts);
%     pattern = {[1, 0, gamma_level, gamma_level], [gamma_level, 0, 1, loop_level],  [0,gamma_level, loop_level, 0] };
% 
%     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
%     %[obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, gamma_dim, 1);
%     %[obj, ln_prefact] = solve_non_lin_and_assign(obj, map, pattern, ln_prefact, nopts);
% 
%     %left up
%     [map, ~] = create_map([0, 1, 0, 0;
%                         0, 2, 3, 6;
%                         0, 4, 5, 0; ], obj.numopts);
%     pattern = {[0, 1, gamma_level, gamma_level]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 1, 0, 0;
%                         7, 2, 3, 6;
%                         0, 4, 5, 0; ], obj.numopts);
%     pattern = {[1, 1, gamma_level, gamma_level]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
% 
%     %up right
%     [map, ~] = create_map([0, 0, 6, 0;
%                         1, 2, 3, 0;
%                         0, 4, 5, 0; ], obj.numopts);
%     pattern = {[gamma_level, 1, 0, loop_level]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 0, 6, 0;
%                         1, 2, 3, 7;
%                         0, 4, 5, 0; ], obj.numopts);
%     pattern = {[gamma_level, 1, 0, loop_level]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
% 
%     %%%%%%%%%%
%     %%%%%%right
%     [map, ~] = create_map([0, 1, 0;
%                         2, 3, 0;
%                         4, 5, 0;
%                         0, 6, 0], obj.numopts);
%     pattern = { [gamma_level +1, 1, 0, gamma_level +1], [loop_level, gamma_level +1, 0, 1], [0,0,gamma_level +1,loop_level]};
%     
%     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
%     %[obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %right up
%     [map, ~] = create_map([0, 0, 0;
%                         2, 3, 1;
%                         4, 5, 0;
%                         0, 6, 0], obj.numopts);
%     pattern = {[gamma_level +1, 0, 1, gamma_level +1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 7, 0;
%                         2, 3, 1;
%                         4, 5, 0;
%                         0, 6, 0], obj.numopts);
%     pattern = {[gamma_level +1, 1, 1, gamma_level +1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %right down
%     [map, ~] = create_map([0, 1, 0;
%                         2, 3, 0;
%                         4, 5, 6;
%                         0, 0, 0], obj.numopts);
%     pattern = {[loop_level, gamma_level +1, 1, 0]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 1, 0;
%                         2, 3, 0;
%                         4, 5, 6;
%                         0, 7, 0], obj.numopts);
%     pattern = {[loop_level, gamma_level +1, 1, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
% 
%     %%%%%%%%%%
%     %%%%%%down
%     [map, ~] = create_map([0, 2, 3, 0;
%                         1, 4, 5, 6; ], obj.numopts);
%     pattern = {[1, beta_level + 1, gamma_level, 0], [gamma_level, beta_level, 1, 0]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %down left
%     [map, ~] = create_map([0, 2, 3, 0;
%                         0, 4, 5, 6;
%                         0, 1, 0, 0], obj.numopts);
%     pattern = {[0, beta_level + 1, gamma_level, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 2, 3, 0;
%                         7, 4, 5, 6;
%                         0, 1, 0, 0], obj.numopts);
%     pattern = {[1, beta_level + 1, gamma_level, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %down right
%     [map, ~] = create_map([0, 2, 3, 0;
%                         1, 4, 5, 0;
%                         0, 0, 6, 0], obj.numopts);
%     pattern = {[gamma_level, beta_level, 0, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 2, 3, 0;
%                         1, 4, 5, 7;
%                         0, 0, 6, 0], obj.numopts);
%     pattern = {[gamma_level, beta_level, 1, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %%%%%%%%%%
%     %%%%%%left
%     [map, ~] = create_map([0, 1, 0;
%                         0, 2, 3;
%                         0, 4, 5;
%                         0, 6, 0], obj.numopts);
%     pattern = {[0, 1, beta_level, gamma_level], [0, gamma_level, beta_level + 1, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %left up
%     [map, ~] = create_map([0, 0, 0;
%                         1, 2, 3;
%                         0, 4, 5;
%                         0, 6, 0], obj.numopts);
%     pattern = {[1, 0, beta_level, gamma_level]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 7, 0;
%                         1, 2, 3;
%                         0, 4, 5;
%                         0, 6, 0], obj.numopts);
%     pattern = {[1, 1, beta_level, gamma_level]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %left down
%     [map, ~] = create_map([0, 1, 0;
%                         0, 2, 3;
%                         6, 4, 5;
%                         0, 0, 0], obj.numopts);
%     pattern = {[1, gamma_level, beta_level + 1, 0]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
%     %
%     [map, ~] = create_map([0, 1, 0;
%                         0, 2, 3;
%                         6, 4, 5;
%                         0, 7, 0], obj.numopts);
%     pattern = {[1, gamma_level, beta_level + 1, 1]};
%     [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact, -1, 1);
% 
%     if testing == 1
%         calculate_error(obj, [1, 2, 3, 6; 0, 4, 5, 0; ], obj.numopts)
%         calculate_error(obj, [0, 1, 0, 0;0, 2, 3, 6;0, 4, 5, 0; ], obj.numopts)
%         calculate_error(obj, [0, 1, 0, 0;7, 2, 3, 6;0, 4, 5, 0; ], obj.numopts)
%     end

    %%
end
