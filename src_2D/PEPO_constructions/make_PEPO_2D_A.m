function obj = makePEPO3(obj)
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

    part = obj.get_middle_part(...
        {[], [], [], []}, [1, 2], 0);

    [U, S, V] = svd(reshape(part, d^2, d^2));

    sqrt_S = diag(diag(S).^0.5);

    block_l = permute(reshape(U * sqrt_S, [1, d, d, d^2]), [2, 3, 1, 4]);
    block_r = permute(reshape(sqrt_S * V', [d^2, d, d, 1]), [2, 3, 1, 4]);

    obj.PEPO_cell{1, 1, 2, 1} = reshape(block_l, [d, d, 1, 1, d^2, 1]); %right
    obj.PEPO_cell{1, 1, 1, 2} = reshape(block_l, [d, d, 1, 1, 1, d^2]); %down

    obj.PEPO_cell{2, 1, 1, 1} = reshape(block_r, [d, d, d^2, 1, 1, 1]); %left
    obj.PEPO_cell{1, 2, 1, 1} = reshape(block_r, [d, d, 1, d^2, 1, 1]); %up

    %%%%%%%%%%%%%%%%%create 0--|--1--|--1--|--0 and variants
    if obj.max_index >= 1

        obj.current_max_index = 1;

        obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^4]
        obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^4]

        %solve with new solver;
        map = create_map([1, 2, 3], obj.numopts);

        Tensor = obj.H_exp(map, obj.nf) - ...
            obj.contract_network(map, struct('max_index', obj.current_max_index));

        Tensor_site = reshape(permute(Tensor, site_ordering_permute(map.N)), ...
            [d^2, d^2, d^2]);

        opts.tol = 1e-13;
        opts.maxit = 1;
        opts.print_level = 1;
        opts.optim = [2]; %only optimize 2;
        opts.get_elem_num = [1, 2, 3];
        opts.solve_type = {'', 'matrix_inv', ''};

        elem_list = cell(3, 1);
        elem_list{1} = obj.PEPO_cell{1, 1, 2, 1};
        elem_list{2} = rand(d, d, d^2, 1, d^2, 1);
        elem_list{3} = obj.PEPO_cell{2, 1, 1, 1};

        [elems, err, A] = tensor_ring(elem_list, map, Tensor_site, opts);

        %
        block_11 = elems{2};
        %block_11_bis = obj.get_middle_part(
        %{[1,2],[],[],[2;3]},[1,2;0,3] ); same

        %6E-8
        %             obj.calculate_error( create_map([1 2 3],1)) %  0
        %             obj.calculate_error( create_map([1 2;0 3],1)) % 0
        %             obj.calculate_error( create_map([1 0; 2 3],1)) %2E-16
        %             obj.calculate_error( create_map([1; 2 ;3;],1)) %  0

        %al same block because up and left blocks are the same
        obj.PEPO_cell{2, 1, 2, 1} = reshape(block_11, [d, d, d^2, 1, d^2, 1]);
        obj.PEPO_cell{2, 1, 1, 2} = reshape(block_11, [d, d, d^2, 1, 1, d^2]);
        obj.PEPO_cell{1, 2, 2, 1} = reshape(block_11, [d, d, 1, d^2, d^2, 1]);
        obj.PEPO_cell{1, 2, 1, 2} = reshape(block_11, [d, d, 1, d^2, 1, d^2]);
        %
        %

        if obj.testing == 1
            obj.calculate_error(create_map([1 2 3], obj.numopts))
            obj.calculate_error(create_map([1 2; 0 3], obj.numopts))
            obj.calculate_error(create_map([1 0; 2 3], obj.numopts))
            obj.calculate_error(create_map([1; 2; 3; ], obj.numopts))
        end

        %make cyclic properties better with null space vectors
        opts_h_cyclic_numbered.numbered = 1;
        opts_h_cyclic_numbered.v_cyclic = 0;
        opts_h_cyclic_numbered.h_cyclic = 1;

        map_c = create_map([1, 2, 3], opts_h_cyclic_numbered); %cyclic ring

        Tensor_c = obj.H_exp(map_c, obj.nf) - ...
            obj.contract_network(map_c, struct('max_index', obj.current_max_index));

        Tensor_site_c = reshape(permute(Tensor_c, site_ordering_permute(map.N)), ...
            [d^2, d^2, d^2]);

        opts.tol = 1e-13;
        opts.maxit = 1;
        opts.print_level = 1;
        opts.get_elem_num = [1; 1; 1];
        opts.solve_type = {'fsolve', '', ''};
        opts.null_space

        elem_list{1} = rand(d, d, d^2, 1, d^2, 1);
        [elems, err, A] = tensor_ring(elem_list, map, Tensor_site, opts);

        obj.PEPO_cell{3, 1, 3, 1} = elems{1};

        if obj.testing == 1
            obj.calculate_error(create_map([1 2 3], obj.numopts))%no improvement
            obj.calculate_error(create_map([1 2 3], opts_h_cyclic_numbered))
            obj.calculate_error(create_map([1 2 3 4 5 6], opts_h_cyclic_numbered))
        end

    end

    %%%%%%%%%%%

    obj = obj.cell2matrix(); %save matrix form

end
