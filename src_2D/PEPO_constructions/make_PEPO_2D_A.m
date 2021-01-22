function obj = make_PEPO_2D_A(obj)
    d = obj.dim;

    %todo do this in code
    obj.virtual_level_sizes_horiz = [d^0, d^2, d^2];
    obj.virtual_level_sizes_vert = [d^0, d^2, d^2];

    %%%%%%%%%%single site
    O_0000 = eye(d); %expm( 0*obj.H_1_tensor );
    %obj.nf = trace(O_0000);

    obj.PEPO_cell{1, 1, 1, 1} = reshape(O_0000 / exp(obj.nf), [d, d, 1, 1, 1, 1]);


    %%%%%%%%%%%%%% 0--|--1--|--0 and all other veriants
    obj.current_max_index = 0;

    obj.PEPO_cell{1, 1, 2, 1} = reshape(block_l, [d, d, 1, 1, d^2, 1]); %right
    obj.PEPO_cell{1, 1, 1, 2} = reshape(block_l, [d, d, 1, 1, 1, d^2]); %down

    obj.PEPO_cell{2, 1, 1, 1} = reshape(block_r, [d, d, d^2, 1, 1, 1]); %left
    obj.PEPO_cell{1, 2, 1, 1} = reshape(block_r, [d, d, 1, d^2, 1, 1]); %up

    %al same block because up and left blocks are the same
    obj.PEPO_cell{2, 1, 2, 1} = reshape(block_11, [d, d, d^2, 1, d^2, 1]);
    obj.PEPO_cell{2, 1, 1, 2} = reshape(block_11, [d, d, d^2, 1, 1, d^2]);
    obj.PEPO_cell{1, 2, 2, 1} = reshape(block_11, [d, d, 1, d^2, d^2, 1]);
    obj.PEPO_cell{1, 2, 1, 2} = reshape(block_11, [d, d, 1, d^2, 1, d^2]);


    if obj.testing == 1
        obj.calculate_error(create_map([1 2 3], obj.numopts))
        obj.calculate_error(create_map([1 2; 0 3], obj.numopts))
        obj.calculate_error(create_map([1 0; 2 3], obj.numopts))
        obj.calculate_error(create_map([1; 2; 3; ], obj.numopts))
    end


    obj.PEPO_cell{3, 1, 3, 1} = 0;%elems{1};

    if obj.testing == 1
        obj.calculate_error(create_map([1 2 3], obj.numopts))%no improvement
        obj.calculate_error(create_map([1 2 3], opts_h_cyclic_numbered))
        obj.calculate_error(create_map([1 2 3 4 5 6], opts_h_cyclic_numbered))
    end
end
