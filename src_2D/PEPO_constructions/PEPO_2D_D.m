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

    %% 1 2 alpha_level and 3 1 levels
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

    %% 2 2 levels + 1 alpha_level
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

    %% 2 2 levels + 2 1 alpha_level
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

    %% 3 2 legs and 1 leg %works but takes a few second
    [map, ~] = create_map([0, 0, 8, 0, 0;
                        1, 2, 3, 4, 5;
                        0, 0, 6, 0, 0
                        0, 0, 7, 0, 0], obj.numopts);
    pattern = {[2, 1, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 8, 0, 0;
                        0, 0, 3, 0, 0;
                        1, 2, 6, 4, 5
                        0, 0, 7, 0, 0], obj.numopts);
    pattern = {[2, 2, 2, 1]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 1, 0, 0;
                        0, 2, 0, 0;
                        8, 3, 4, 7;
                        0, 5, 0, 0;
                        0, 6, 0, 0], obj.numopts);
    pattern = {[1, 2, 2, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    [map, ~] = create_map([0, 0, 1, 0;
                        0, 0, 2, 0;
                        7, 3, 4, 8;
                        0, 0, 5, 0;
                        0, 0, 6, 0], obj.numopts);
    pattern = {[2, 2, 1, 2]};
    [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);

    if obj.testing == 1
        calculate_error(obj, [0, 0, 8, 0, 0; 1, 2, 3, 4, 5; 0, 0, 6, 0, 0; 0, 0, 7, 0, 0], obj.numopts)
    end
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

    fprintf("\n")

    obj = cell2matrix(obj);

    alpha_level = obj.current_max_index + 1;
    beta_level = alpha_level + 1;

    obj.current_max_index = obj.current_max_index + 3;
    obj.max_index = obj.current_max_index;
    alpha_dim = 10;
    beta_dim = 10;

    lnlopts = struct('Display', 1, 'maxit', 5);
    %obj.cycle_index = alpha_level;

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, alpha_dim,alpha_dim, beta_dim, beta_dim];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, alpha_dim,alpha_dim, beta_dim, beta_dim];

    %simple loop

    [alpha_map, ~] = create_map([1, 2; 3, 4], obj.numopts);
    alpha_pattern = {[alpha_level, alpha_level, 0, 0], [alpha_level, 0, 0, alpha_level], [0, alpha_level, alpha_level, 0], [0, 0, alpha_level, alpha_level]};
    obj = solve_lin_non_lin_and_assign(obj, alpha_map, alpha_pattern, ln_prefact, lnlopts);
    obj = rescale_PEPO_pattern(obj, alpha_pattern);

    if obj.testing == 1
        calculate_error(obj, [1, 2; 3, 4], obj.numopts)
    end

    % ofset double loop
    onceopts = struct('Display', 1, 'maxit', 2);
    repopts = struct('Display', 1, 'maxit', 20);

    [map1, ~] = create_map([0, 2, 1;
                        6, 4, 3;
                        5, 7, 0], obj.numopts);

    %neglect new 1 loop, will be fixed later

    %
    %     pattern1 = {[beta_level, beta_level, beta_level + 1, beta_level + 1], ...
    %                 [0, 0, beta_level, beta_level], ...
    %                 [beta_level + 1, beta_level+1, 0, 0], ...
    %                 [0, 0, beta_level, beta_level], ...
    %                 [beta_level+1, beta_level + 1, 0, 0], ...
    %                 [0, beta_level, beta_level+1, 0], ...
    %                 [ beta_level,0,0, beta_level+1], ...
    %                 };
    %
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts, pattern1);
    %     obj = rescale_PEPO_pattern(obj, pattern1);
    %
    %     redo simple loop
    %     [alpha_map, ~] = create_map([1, 2; 3, 4], obj.numopts);
    %     alpha_pattern = {[alpha_level, alpha_level, 0, 0], [alpha_level, 0, 0, alpha_level], [0, alpha_level, alpha_level, 0], [0, 0, alpha_level, alpha_level]};
    %     obj = solve_lin_non_lin_and_assign(obj, alpha_map, alpha_pattern, ln_prefact, lnlopts);
    %     obj = rescale_PEPO_pattern(obj, alpha_pattern);
    %
    %     do other loop
    %      [map1, ~] = create_map([1, 2, 0;
    %                         3, 4, 6;
    %                         0, 7, 5; ], obj.numopts);
    %
    %     pattern1 = {[0, beta_level, beta_level, 0], ...
    %                 [beta_level, 0,0, beta_level], ...
    %                 [0, beta_level + 1, beta_level+1, 0], ...
    %                 [beta_level+1, 0, 0, beta_level+1], ...
    %                 };
    %
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %     obj = rescale_PEPO_pattern(obj, pattern1);
    %
    %
    %     if obj.testing == 1
    %         calculate_error(obj, [0, 2, 1;
    %                         6, 4, 3;
    %                         5, 7, 0], obj.numopts)
    %         calculate_error(obj, [1, 2, 0;
    %                         3, 4, 6;
    %                         0, 7, 5; ], obj.numopts)
    %     end
    %
    %     %double loop
    %     [map1, ~] = create_map([1, 2, 5; 3, 4, 6], obj.numopts);
    %     pattern1 = {[beta_level, beta_level, beta_level + 1, 0], [beta_level, 0, beta_level + 1, beta_level]};
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %     obj = rescale_PEPO_pattern(obj, pattern1);
    %     obj = scale_PEPO_pattern(obj, pattern, 0.5)%only solve half of the prob
    %
    %     pattern1 = {[beta_level, beta_level + 1, beta_level + 1, 0], [beta_level, 0, beta_level + 1, beta_level + 1]};
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %     obj = rescale_PEPO_pattern(obj, pattern1);
    %
    %     pattern1 = {[beta_level, beta_level, beta_level + 1, 0], [beta_level, 0, beta_level + 1, beta_level], ...
    %                 [beta_level, beta_level + 1, beta_level + 1, 0], [beta_level, 0, beta_level + 1, beta_level + 1]};
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %
    %     obj = scale_PEPO_pattern(obj, pattern1, 1/2);
    %     double loop
    %     [map1, ~] = create_map([1, 2; 3, 4; 5, 6], obj.numopts);
    %
    %     pattern1 = {[0, beta_level, beta_level, beta_level + 1], [beta_level, beta_level, 0, beta_level + 1]}
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, onceopts);
    %     obj = rescale_PEPO_pattern(obj, pattern1);
    %     obj = scale_PEPO_pattern(obj, pattern, 0.5)%only solve half of the prob
    %
    %     pattern1 = {[0, beta_level, beta_level + 1, beta_level + 1], [beta_level + 1, beta_level, 0, beta_level + 1]};
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, onceopts);
    %     obj = rescale_PEPO_pattern(obj, pattern1);
    %
    %     pattern1 = {[0, beta_level, beta_level + 1, beta_level + 1], [beta_level + 1, beta_level, 0, beta_level + 1], ...
    %                 [0, beta_level, beta_level + 1, beta_level + 1], [beta_level + 1, beta_level, 0, beta_level + 1]};
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %
    %     if obj.testing == 1
    %         calculate_error(obj, [1, 2, 5; 3, 4, 6], obj.numopts)
    %         calculate_error(obj, [1, 2; 3, 4; 5, 6], obj.numopts)
    %     end

    %     % ofset double loop
    %     onceopts = struct('Display', 1, 'maxit', 2);
    %     repopts = struct('Display', 1, 'maxit', 10);
    %
    %     [map1, ~] = create_map([1, 2, 0;
    %                         3, 4, 6;
    %                         0, 7, 5; ], obj.numopts);
    %
    %     pattern1 = {[beta_level, beta_level, beta_level + 1, beta_level + 1], ...
    %                 [beta_level + 1, 0, 0, alpha_level], ...
    %                 [0, alpha_level, beta_level, 0], ...
    %                 [0, beta_level + 1, alpha_level, 0], ...
    %                 [alpha_level, 0, 0, beta_level], ...
    %                 };
    %
    %     %obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %     %obj = rescale_PEPO_pattern(obj, pattern1);
    %
    %     [map2, ~] = create_map([0, 2, 1;
    %                         6, 4, 3;
    %                         5, 7, 0], obj.numopts);
    %
    %     pattern2 = {[0, 0, beta_level, alpha_level], ...
    %                 [beta_level + 1, alpha_level, 0, 0], ...
    %                 [0, 0, alpha_level, beta_level], ...
    %                 [alpha_level, beta_level + 1, 0, 0], ...
    %                 };
    %
    %     map = {map1,map2};
    %     pattern = [pattern1,pattern2];
    %
    %
    %     [obj, ln_prefact] = solve_non_lin_and_assign(obj, map, pattern, ln_prefact, struct('Display', 'iter-detailed'))     ;
    %
    %     %obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %     %obj = rescale_PEPO_pattern(obj, pattern1);
    %
    %
    %     if obj.testing == 1
    %         calculate_error(obj, [0, 2, 1;
    %                         6, 4, 3;
    %                         5, 7, 0], obj.numopts)
    %         calculate_error(obj, [1, 2, 0;
    %                         3, 4, 6;
    %                         0, 7, 5; ], obj.numopts)
    %     end
    %
    %     % %double loop
    %     [map1, ~] = create_map([1, 2, 5; 3, 4, 6], obj.numopts);
    %     pattern1 = {[beta_level, beta_level, beta_level + 1, 0], [beta_level, 0, beta_level + 1, beta_level]};
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %     obj = rescale_PEPO_pattern(obj, pattern1);
    %     obj = scale_PEPO_pattern(obj, pattern, 0.5)%only solve half of the prob
    %
    %     pattern1 = {[beta_level, beta_level + 1, beta_level + 1, 0], [beta_level, 0, beta_level + 1, beta_level + 1]};
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %     obj = rescale_PEPO_pattern(obj, pattern1);
    %
    %     pattern1 = {[beta_level, beta_level, beta_level + 1, 0], [beta_level, 0, beta_level + 1, beta_level], ...
    %                 [beta_level, beta_level + 1, beta_level + 1, 0], [beta_level, 0, beta_level + 1, beta_level + 1]};
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %
    %     %obj = scale_PEPO_pattern(obj, pattern1, 1/2);
    %     %double loop
    %     [map1, ~] = create_map([1, 2; 3, 4; 5, 6], obj.numopts);
    %
    %     pattern1 = {[0, beta_level, beta_level, beta_level + 1], [beta_level, beta_level, 0, beta_level + 1]}
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, onceopts);
    %     obj = rescale_PEPO_pattern(obj, pattern1);
    %     obj = scale_PEPO_pattern(obj, pattern, 0.5)%only solve half of the prob
    %
    %     pattern1 = {[0, beta_level, beta_level + 1, beta_level + 1], [beta_level + 1, beta_level, 0, beta_level + 1]};
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, onceopts);
    %     obj = rescale_PEPO_pattern(obj, pattern1);
    %
    %     pattern1 = {[0, beta_level, beta_level + 1, beta_level + 1], [beta_level + 1, beta_level, 0, beta_level + 1], ...
    %                 [0, beta_level, beta_level + 1, beta_level + 1], [beta_level + 1, beta_level, 0, beta_level + 1]};
    %     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, repopts);
    %
    %     if obj.testing == 1
    %         calculate_error(obj, [1, 2, 5; 3, 4, 6], obj.numopts)
    %         calculate_error(obj, [1, 2; 3, 4; 5, 6], obj.numopts)
    %     end

    % % %double loop
    % [map1, ~] = create_map([1, 2, 5; 3, 4, 6], obj.numopts);

    % pattern1 = {[beta_level, beta_level, beta_level + 1, 0], [beta_level, 0, beta_level + 1, beta_level]};
    % pattern2 = {[beta_level, beta_level + 1, beta_level + 1, 0], [beta_level, 0, beta_level + 1, beta_level + 1]};
    % corners = {[beta_level + 1, 0, 0, beta_level + 1], [beta_level + 1, beta_level + 1, 0, 0], ...
    %             [0, beta_level, beta_level, 0], [0, 0, beta_level, beta_level]};

    % pattern = [pattern1, corners];
    % obj = solve_lin_non_lin_and_assign(obj, map1, pattern, ln_prefact, onceopts);
    % obj = rescale_PEPO_pattern(obj, pattern);
    % obj = scale_PEPO_pattern(obj, pattern, 0.5);

    % obj = solve_lin_non_lin_and_assign(obj, map1, pattern2, ln_prefact, onceopts);
    % obj = rescale_PEPO_pattern(obj, pattern2);

    % pattern = [pattern1, pattern2];
    % obj = solve_lin_non_lin_and_assign(obj, map1, pattern, ln_prefact, repopts);

    % %obj = scale_PEPO_pattern(obj, pattern1, 1/2);
    % %double loop
    % [map1, ~] = create_map([1, 2; 3, 4; 5, 6], obj.numopts);

    % pattern1 = {[0, beta_level, beta_level, beta_level + 1], [beta_level, beta_level, 0, beta_level + 1]};
    % pattern2 = {[0, beta_level, beta_level + 1, beta_level + 1], [beta_level + 1, beta_level, 0, beta_level + 1]};
    % corners = {[0, beta_level + 1, beta_level + 1, 0], [beta_level, 0, 0, beta_level]};

    % pattern = [pattern1, corners];
    % obj = solve_lin_non_lin_and_assign(obj, map1, pattern, ln_prefact, onceopts);
    % obj = rescale_PEPO_pattern(obj, pattern);
    % obj = scale_PEPO_pattern(obj, pattern, 0.5);

    % obj = solve_lin_non_lin_and_assign(obj, map1, pattern2, ln_prefact, onceopts);
    % obj = rescale_PEPO_pattern(obj, pattern2);

    % pattern = [pattern1, pattern2];
    % obj = solve_lin_non_lin_and_assign(obj, map1, pattern, ln_prefact, repopts);

    % if obj.testing == 1
    %     calculate_error(obj, [1, 2, 5; 3, 4, 6], obj.numopts)
    %     calculate_error(obj, [1, 2; 3, 4; 5, 6], obj.numopts)
    % end

    % %double loop
    [map1, ~] = create_map([1, 2, 5; 3, 4, 6], obj.numopts);

    pattern1 = {[beta_level, alpha_level, beta_level + 1, 0], [beta_level, 0, beta_level + 1, alpha_level]};
    corners = {[beta_level + 1, 0, 0, beta_level + 1], [beta_level + 1, beta_level + 1, 0, 0], ...
                [0, beta_level, beta_level, 0], [0, 0, beta_level, beta_level]};

    pattern = [pattern1, corners];
    obj = solve_lin_non_lin_and_assign(obj, map1, pattern, ln_prefact, repopts);
    obj = rescale_PEPO_pattern(obj, pattern);

    %obj = scale_PEPO_pattern(obj, pattern1, 1/2);
    %double loop
    [map1, ~] = create_map([1, 2; 3, 4; 5, 6], obj.numopts);

    pattern1 = {[0, beta_level, alpha_level, beta_level + 1], [alpha_level, beta_level, 0, beta_level + 1]};
    corners = {[0, beta_level + 1, beta_level + 1, 0], [beta_level, 0, 0, beta_level]};

    pattern = [pattern1, corners];
    obj = solve_lin_non_lin_and_assign(obj, map1, pattern, ln_prefact, repopts);
    obj = rescale_PEPO_pattern(obj, pattern);

    if obj.testing == 1
        calculate_error(obj, [1, 2, 5; 3, 4, 6], obj.numopts)
        calculate_error(obj, [1, 2; 3, 4; 5, 6], obj.numopts)
    end

    % ofset double loop
    onceopts = struct('Display', 1, 'maxit', 1);
    repopts = struct('Display', 1, 'maxit', 10);

    % %
%     [map1, ~] = create_map([1, 2, 0;
%                         3, 4, 6;
%                         0, 7, 5; ], obj.numopts);
% 
%     pattern1 = {[alpha_level, beta_level, alpha_level, beta_level + 1], ...
%                 [alpha_level, 0, 0, beta_level + 1], ...
%                 [0, beta_level, alpha_level, 0]};
% 
%     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, onceopts);
%     obj = rescale_PEPO_pattern(obj, pattern1);
%     %
%     [map1, ~] = create_map([0, 2, 1;
%                         6, 4, 3;
%                         5, 7, 0], obj.numopts);
% 
%     pattern1 = {[beta_level, alpha_level, beta_level + 1, alpha_level], ...
%                 };
%     %[0, 0, beta_level + 1, alpha_level], ...
%     %[beta_level, alpha_level, 0, 0]};
% 
%     obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, onceopts);
%     obj = rescale_PEPO_pattern(obj, pattern1);

    % %
    % [map1, ~] = create_map([1, 2, 0;
    %                     3, 4, 6;
    %                     0, 7, 5; ], obj.numopts);

    % pattern1 = {[alpha_level, beta_level, beta_level + 1, beta_level + 1], ...
    %             [0, beta_level, alpha_level, 0]};

    % obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, onceopts);
    % obj = rescale_PEPO_pattern(obj, pattern1);
    % %
    % [map1, ~] = create_map([0, 2, 1;
    %                     6, 4, 3;
    %                     5, 7, 0], obj.numopts);

    % pattern1 = {[beta_level, beta_level, beta_level + 1, alpha_level]};
    % corners = {...
    %             [beta_level, 0, 0, alpha_level], ...
    %             [beta_level + 1, alpha_level, 0, 0], ...
    %             };

    % %[0, 0, beta_level + 1, alpha_level], ...
    % %[beta_level, alpha_level, 0, 0]};

    % obj = fill_rand(obj, corners, mean(obj.PEPO_cell{1, beta_level + 1, alpha_level + 1, 1}, 'all'));

    % %redo simple loop
    % [alpha_map, ~] = create_map([1, 2; 3, 4], obj.numopts);
    % alpha_pattern = {[alpha_level, alpha_level, 0, 0], [alpha_level, 0, 0, alpha_level], [0, alpha_level, alpha_level, 0], [0, 0, alpha_level, alpha_level]};
    % [obj,~] = solve_non_lin_and_assign(obj, {alpha_map}, alpha_pattern, ln_prefact, struct(struct('Display','iter-detailed')));

    % %obj = solve_lin_non_lin_and_assign(obj, alpha_map, alpha_pattern, ln_prefact, lnlopts);
    % %obj = rescale_PEPO_pattern(obj, alpha_pattern);

    % %continue
    % obj = solve_lin_and_assign(obj, map1, pattern1, ln_prefact, -1, 1);
    % %obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, onceopts);
    % %obj = rescale_PEPO_pattern(obj, pattern1);

    
            %
%                       %
%         [map1, ~] = create_map([0, 2, 1;
%                             6, 4, 3;
%                             5, 7, 0], obj.numopts);
%     
%         pattern1 = {[alpha_level, alpha_level+1, alpha_level+1, alpha_level ], ...
%                     [0,0,alpha_level+1, alpha_level+1], ...
%                     [alpha_level+1, 0, 0, alpha_level+1], ...
%                     [alpha_level+1, alpha_level+1 0, 0,], ...
%                     };
%                 
%         obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, onceopts);
%         obj = rescale_PEPO_pattern(obj, pattern1);
%             
%            % 
%         [map1, ~] = create_map([1, 2, 0;
%                     3, 4, 6;
%                     0, 7, 5; ], obj.numopts);
%     
%         pattern1 = {[beta_level, beta_level, beta_level + 1, beta_level + 1]};
%     
%         obj = solve_lin_non_lin_and_assign(obj, map1, pattern1, ln_prefact, onceopts);
%         obj = rescale_PEPO_pattern(obj, pattern1);

  

        %
    
    
        %[0, 0, beta_level + 1, alpha_level], ...
        %[beta_level, alpha_level, 0, 0]};
    
    if obj.testing == 1
        calculate_error(obj, [1, 2, 0;
                        3, 4, 6;
                        0, 7, 5; ], obj.numopts)

        calculate_error(obj, [0, 2, 1;
                        6, 4, 3;
                        5, 7, 0], obj.numopts)
    end

end
