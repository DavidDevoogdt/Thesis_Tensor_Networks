function obj = make_PEPO_2D_A(obj)
    d = obj.dim;
    ln_prefact = obj.nf;

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

    %% 1 2 a and 3 1 levels
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

    %% 2 2 levels + 1 a
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

    %% 2 2 levels + 2 1 a
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

    %         %% 3 2 legs and 1 leg %works but takes a few second
    %         [map, ~] = create_map([0, 0, 8, 0, 0;
    %                             1, 2, 3, 4, 5;
    %                             0, 0, 6, 0, 0
    %                             0, 0, 7, 0, 0], obj.numopts);
    %         pattern = {[2, 1, 2, 2]};
    %         [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %         %
    %         [map, ~] = create_map([0, 0, 8, 0, 0;
    %                             0, 0, 3, 0, 0;
    %                             1, 2, 6, 4, 5
    %                             0, 0, 7, 0, 0], obj.numopts);
    %         pattern = {[2, 2, 2, 1]};
    %         [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %         %
    %         [map, ~] = create_map([0, 1, 0, 0;
    %                             0, 2, 0, 0;
    %                             8, 3, 4, 7;
    %                             0, 5, 0, 0;
    %                             0, 6, 0, 0], obj.numopts);
    %         pattern = {[1, 2, 2, 2]};
    %         [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %         %
    %         [map, ~] = create_map([0, 0, 1, 0;
    %                             0, 0, 2, 0;
    %                             7, 3, 4, 8;
    %                             0, 0, 5, 0;
    %                             0, 0, 6, 0], obj.numopts);
    %         pattern = {[2, 2, 1, 2]};
    %         [obj, ~, ~, ln_prefact, ~] = solve_lin_and_assign(obj, map, pattern, ln_prefact);
    %
    %         if obj.testing == 1
    %             calculate_error(obj, [0, 0, 8, 0, 0; 1, 2, 3, 4, 5; 0, 0, 6, 0, 0; 0, 0, 7, 0, 0], obj.numopts)
    %         end
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

    %fprintf("\n")

    obj = cell2matrix(obj);

    %     x = 3;
    %     y = x + 1;
    %     a = y + 1;
    %     b = a + 1;
    %
    %     %obj.PEPO_cell{zero_level+1,zero_level+1,zero_level+1,zero_level+1}= obj.PEPO_cell{1,1,1,1};
    %
    %     alpha_dim = 19;
    %     beta_dim = 19;
    %
    %     lnlopts = struct('Display', 1, 'maxit', 1);
    %     %obj.cycle_index = a;
    %
    %     obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^2, d^2 alpha_dim, beta_dim];
    %     obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^2, d^2, alpha_dim, beta_dim];
    %
    %  f = 1;
    %     obj.PEPO_cell{y + 1, 1, 1, 1} = reshape(f * eye(obj.dim^2) / exp(obj.nf), [d, d, d^2, 1, 1, 1]);
    %     obj.PEPO_cell{1, y + 1, 1, 1} = reshape(f * eye(obj.dim^2) / exp(obj.nf), [d, d, 1, d^2, 1, 1]);
    %     obj.PEPO_cell{1, 1, x + 1, 1} = reshape(f * eye(obj.dim^2) / exp(obj.nf), [d, d, 1, 1, d^2, 1]);
    %     obj.PEPO_cell{1, 1, 1, x + 1} = reshape(f * eye(obj.dim^2) / exp(obj.nf), [d, d, 1, 1, 1, d^2]);
    %

    x = 1;
    y = 1;
    a = 3;
    b = a + 1;

    %obj.PEPO_cell{zero_level+1,zero_level+1,zero_level+1,zero_level+1}= obj.PEPO_cell{1,1,1,1};

    alpha_dim = 10;
    beta_dim = 16;

    lnlopts = struct('Display', 0, 'maxit', 1);
    %obj.cycle_index = a;

    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, alpha_dim, beta_dim];
    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, alpha_dim, beta_dim];

    %obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^2, d^2 alpha_dim];
    %obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^2, d^2, alpha_dim];

    obj.current_max_index = numel(obj.virtual_level_sizes_horiz);
    obj.max_index = obj.current_max_index;

    %simple loop

    [alpha_map, ~] = create_map([1, 2; 3, 4], obj.numopts);
    alpha_pattern = {[a, a, 0, 0], ...
                    [a, 0, 0, a], ...
                    [0, a, a, 0], ...
                    [0, 0, a, a]};

    obj = solve_lin_non_lin_and_assign(obj, alpha_map, alpha_pattern, ln_prefact, lnlopts);
    obj = rescale_PEPO_pattern(obj, alpha_pattern);

    if obj.testing == 1
        calculate_error(obj, [1, 2; 3, 4], obj.numopts)
    end

    %% double loops
    %[x,x,y,y]

    dbool = 0;

    [map, ~] = create_map([
                        0, 1, 0, 0;
                        0, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 0, 0, 0]);
    pattern = {[0, x, a, a]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    [map, ~] = create_map([
                        0, 0, 0, 0;
                        1, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 0, 0, 0]);
    pattern = {[x, 0, a, a]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    if dbool == 1
        [map, ~] = create_map([
                            0, 1, 0, 0;
                            1, 1, 1, 0;
                            0, 1, 1, 0;
                            0, 0, 0, 0]);
        pattern = {[x, x, a, a]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    end

    [map, ~] = create_map([
                        0, 0, 1, 0;
                        0, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 0, 0, 0]);
    pattern = {[a, x, 0, a]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    [map, ~] = create_map([
                        0, 0, 0, 0;
                        0, 1, 1, 1;
                        0, 1, 1, 0;
                        0, 0, 0, 0]);
    pattern = {[a, 0, y, a]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    if dbool == 1
        [map, ~] = create_map([
                            0, 0, 1, 0;
                            0, 1, 1, 1;
                            0, 1, 1, 0;
                            0, 0, 0, 0]);
        pattern = {[a, x, y, a]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end
    %
    [map, ~] = create_map([
                        0, 0, 0, 0;
                        0, 1, 1, 0;
                        0, 1, 1, 1;
                        0, 0, 0, 0]);
    pattern = {[a, a, y, 0]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    [map, ~] = create_map([
                        0, 0, 0, 0;
                        0, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 0, 1, 0]);
    pattern = {[a, a, 0, y]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    if dbool == 1
        [map, ~] = create_map([
                            0, 0, 0, 0;
                            0, 1, 1, 0;
                            0, 1, 1, 1;
                            0, 0, 1, 0]);
        pattern = {[a, a, y, y]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end
    %
    [map, ~] = create_map([
                        0, 0, 0, 0;
                        0, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 1, 0, 0]);
    pattern = {[0, a, a, y]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    [map, ~] = create_map([
                        0, 0, 0, 0;
                        0, 1, 1, 0;
                        1, 1, 1, 0;
                        0, 0, 0, 0]);
    pattern = {[x, a, a, 0]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    if dbool == 1
        [map, ~] = create_map([
                            0, 0, 0, 0;
                            0, 1, 1, 0;
                            1, 1, 1, 0;
                            0, 1, 0, 0]);
        pattern = {[x, a, a, y]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end

    if obj.testing == 1
        calculate_error(obj, [0, 0, 0;
                        1, 1, 1;
                        0, 1, 1], struct)

        calculate_error(obj, [0, 1, 0;
                        1, 1, 1;
                        0, 1, 1], struct)
    end

    %% double corrections
    if obj.testing == 1
        fprintf("\n single beta\n")
    end
    %up
    [map, ~] = create_map([
                        0, 0, 1, 0;
                        1, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 0, 0, 0]);
    pattern = {[x, 0, b, a], [b, x, 0, a]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    obj = rescale_PEPO_pattern(obj, pattern);

    [map, ~] = create_map([
                        0, 0, 0, 0;
                        1, 1, 1, 1;
                        0, 1, 1, 0;
                        0, 0, 0, 0]);
    pattern = {[b, 0, y, a]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    if dbool == 1
        [map, ~] = create_map([
                            0, 0, 1, 0;
                            1, 1, 1, 1;
                            0, 1, 1, 0;
                            0, 0, 0, 0]);
        pattern = {[b, x, y, a]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end

    [map, ~] = create_map([
                        0, 1, 0, 0;
                        0, 1, 1, 1;
                        0, 1, 1, 0;
                        0, 0, 0, 0]);
    pattern = {[0, x, b, a]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    if dbool == 1
        [map, ~] = create_map([
                            0, 1, 0, 0;
                            1, 1, 1, 1;
                            0, 1, 1, 0;
                            0, 0, 0, 0]);
        pattern = {[x, x, b, a]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end
    %left%%%%
    [map, ~] = create_map([
                        0, 0, 1, 0;
                        0, 1, 1, 0;
                        0, 1, 1, 1;
                        0, 0, 0, 0]);
    pattern = {[a, x, 0, b], [a, b, y, 0]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    obj = rescale_PEPO_pattern(obj, pattern);

    [map, ~] = create_map([
                        0, 0, 1, 0;
                        0, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 0, 1, 0]);
    pattern = {[a, b, 0, y]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %
    if dbool == 1
        [map, ~] = create_map([
                            0, 0, 1, 0;
                            0, 1, 1, 0;
                            0, 1, 1, 1;
                            0, 0, 1, 0]);
        pattern = {[a, b, y, y]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end

    [map, ~] = create_map([
                        0, 0, 0, 0;
                        0, 1, 1, 1;
                        0, 1, 1, 0;
                        0, 0, 1, 0]);
    pattern = {[a, 0, y, b]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    if dbool == 1
        [map, ~] = create_map([
                            0, 0, 1, 0;
                            0, 1, 1, 1;
                            0, 1, 1, 0;
                            0, 0, 1, 0]);
        pattern = {[a, x, y, b]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end
    %down
    [map, ~] = create_map([
                        0, 0, 0, 0;
                        0, 1, 1, 0;
                        1, 1, 1, 1;
                        0, 0, 0, 0]);
    pattern = {[b, a, y, 0], [x, a, b, 0]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    obj = rescale_PEPO_pattern(obj, pattern);

    [map, ~] = create_map([
                        0, 0, 0, 0;
                        0, 1, 1, 0;
                        0, 1, 1, 1;
                        0, 1, 0, 0]);
    pattern = {[0, a, b, y]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    if dbool == 1

        [map, ~] = create_map([
                            0, 0, 0, 0;
                            0, 1, 1, 0;
                            1, 1, 1, 1;
                            0, 1, 0, 0]);
        pattern = {[x, a, b, y]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end

    [map, ~] = create_map([
                        0, 0, 0, 0;
                        0, 1, 1, 0;
                        1, 1, 1, 0;
                        0, 0, 1, 0]);
    pattern = {[b, a, 0, y]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    if dbool == 1
        [map, ~] = create_map([
                            0, 0, 0, 0;
                            0, 1, 1, 0;
                            1, 1, 1, 1;
                            0, 0, 1, 0]);
        pattern = {[b, a, y, y]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end

    %left
    [map, ~] = create_map([
                        0, 1, 0, 0;
                        0, 1, 1, 0;
                        1, 1, 1, 0;
                        0, 0, 0, 0]);
    pattern = {[x, b, a, 0], [0, x, a, b]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    obj = rescale_PEPO_pattern(obj, pattern);

    [map, ~] = create_map([
                        0, 1, 0, 0;
                        0, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 1, 0, 0]);
    pattern = {[0, b, a, y]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);

    if dbool == 1
        [map, ~] = create_map([
                            0, 1, 0, 0;
                            0, 1, 1, 0;
                            1, 1, 1, 0;
                            0, 1, 0, 0]);
        pattern = {[x, b, a, y]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end

    [map, ~] = create_map([
                        0, 0, 0, 0;
                        1, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 1, 0, 0]);
    pattern = {[x, 0, a, b]};
    obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %
    if dbool == 1
        [map, ~] = create_map([
                            0, 1, 0, 0;
                            1, 1, 1, 0;
                            0, 1, 1, 0;
                            0, 1, 0, 0]);
        pattern = {[x, x, a, b]};
        obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    end

    if obj.testing == 1
        calculate_error(obj, [
                        0, 0, 1, 0;
                        1, 1, 1, 0;
                        0, 1, 1, 0;
                        0, 0, 0, 0], struct)
    end

    %     %%three neighbouring sites
    %     fprintf("\n double beta\n")
    %
    %     dbool = 0;
    %
    %     [map, ~] = create_map([
    %                         0, 0, 1, 0;
    %                         1, 1, 1, 0;
    %                         0, 1, 1, 0;
    %                         0, 0, 1, 0]);
    %     pattern = {[b, x, 0, b]};
    %     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %
    %     [map, ~] = create_map([
    %                         0, 0, 0, 0;
    %                         1, 1, 1, 1;
    %                         0, 1, 1, 0;
    %                         0, 0, 1, 0]);
    %     pattern = {[b, 0, x, b]};
    %     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %
    %     if dbool == 1
    %         [map, ~] = create_map([
    %                             0, 0, 1, 0;
    %                             1, 1, 1, 1;
    %                             0, 1, 1, 0;
    %                             0, 0, 1, 0]);
    %         pattern = {[b, x, x, b]};
    %         obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %     end
    %     %
    %     [map, ~] = create_map([
    %                         0, 0, 1, 0;
    %                         0, 1, 1, 0;
    %                         1, 1, 1, 0;
    %                         0, 0, 1, 0]);
    %     pattern = {[b, b, 0, x]};
    %     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %
    %     [map, ~] = create_map([
    %                         0, 0, 1, 0;
    %                         0, 1, 1, 0;
    %                         1, 1, 1, 1;
    %                         0, 0, 0, 0]);
    %     pattern = {[b, b, x, 0]};
    %     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %
    %     if dbool == 1
    %         [map, ~] = create_map([
    %                             0, 0, 1, 0;
    %                             0, 1, 1, 0;
    %                             1, 1, 1, 1;
    %                             0, 0, 1, 0]);
    %         pattern = {[b, b, x, x]};
    %         obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %     end
    %     %
    %     [map, ~] = create_map([
    %                         0, 1, 0, 0;
    %                         0, 1, 1, 0;
    %                         1, 1, 1, 1;
    %                         0, 0, 0, 0]);
    %     pattern = {[x, b, b, 0]};
    %     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %
    %     [map, ~] = create_map([
    %                         0, 1, 0, 0;
    %                         0, 1, 1, 0;
    %                         0, 1, 1, 1;
    %                         0, 1, 0, 0]);
    %     pattern = {[0, b, b, x]};
    %     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %
    %     if dbool == 1
    %         [map, ~] = create_map([
    %                             0, 1, 0, 0;
    %                             0, 1, 1, 0;
    %                             1, 1, 1, 1;
    %                             0, 1, 0, 0]);
    %         pattern = {[x, b, b, x]};
    %         obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %     end
    %     %
    %     [map, ~] = create_map([
    %                         0, 1, 0, 0;
    %                         0, 1, 1, 1;
    %                         0, 1, 1, 0;
    %                         0, 1, 0, 0]);
    %     pattern = {[0, x, b, b]};
    %     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %
    %     [map, ~] = create_map([
    %                         0, 0, 0, 0;
    %                         1, 1, 1, 1;
    %                         0, 1, 1, 0;
    %                         0, 1, 0, 0]);
    %     pattern = {[x, 0, b, b]};
    %     obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %
    %     if dbool == 1
    %         [map, ~] = create_map([
    %                             0, 1, 0, 0;
    %                             1, 1, 1, 1;
    %                             0, 1, 1, 0;
    %                             0, 1, 0, 0]);
    %         pattern = {[x, x, b, b]};
    %         obj = solve_lin_non_lin_and_assign(obj, map, pattern, ln_prefact, lnlopts);
    %     end
end
