function part = get_middle_part(obj, inv_maps, map, perm)
    %calculates residual hamiltonian and inverts the legs given in
    %inv_maps in the order {left,up,right,down}. left and up
    %should be increasing indices to central pepo cell, right and
    %down decreasing.

    if nargin < 4
        perm = 1;
    end

    d = obj.dim;

    [map, b_map] = create_map(map, obj.numopts);

    Tensor = obj.H_exp(map, obj.nf) - obj.contract_network(map, struct('max_index', obj.current_max_index));

    Tensor_site = permute(Tensor, site_ordering_permute(map.N));

    vector_sizes = [0, 0, 0, 0, 0];

    for l = 1:2
        i_map = inv_maps{l};

        if size(i_map, 1) ~= 0
            vector_sizes(l) = sum(i_map ~= 0) - 1;
        end

    end

    for l = 3:4
        i_map = inv_maps{l};

        if size(i_map, 1) ~= 0
            vector_sizes(l + 1) = sum(i_map ~= 0) - 1;
        end

    end

    vector_sizes(3) = map.N - sum(vector_sizes);
    vector_sizes = d.^(2 * vector_sizes);

    Tensor_site = reshape(Tensor_site, vector_sizes);

    if obj.testing
        Tensor_site_cpy = Tensor_site;
    end

    function [l_chain_site, bond_size] = get_l_chain(l_map, l_map_2)
        l_map_2 = create_map(l_map_2, obj.numopts);
        l_map = create_map(l_map, obj.numopts);

        l_tensors = cell(1, l_map.N - 1);

        last_index = 0;

        for i = 1:l_map.N - 1
            index_set = [0, 0, 0, 0];

            index_set_mask = l_map.leg_list{i}(3:end) > 0;
            non_empty = find(index_set_mask);

            if size(non_empty, 2) == 1
                last_index = non_empty(1);
                index_set(last_index) = i;
            else
                prev_index = mod(last_index + 2, 4);

                if non_empty(1) == prev_index
                    index_set(non_empty(1)) = i - 1;
                    index_set(non_empty(2)) = i;
                    last_index = non_empty(2);
                else
                    index_set(non_empty(2)) = i - 1;
                    index_set(non_empty(1)) = i;
                    last_index = non_empty(1);
                end

            end

            l_tensors{i} = obj.PEPO_cell{index_set(1) + 1, index_set(2) + 1, index_set(3) + 1, index_set(4) + 1};
        end

        leg_list_copy = l_map_2.leg_list;
        l_chain = ncon(l_tensors, leg_list_copy);

        l_chain_size = size(l_chain);
        l_size = d^(2 * l_map_2.N);

        bond_size = numel(l_chain) / l_size;
        l_chain_size_new = [l_chain_size(1:2 * l_map_2.N), bond_size];

        l_chain = reshape(l_chain, l_chain_size_new);

        l_permute = [site_ordering_permute(l_map_2.N); 2 * l_map_2.N + 1];
        l_chain_site = reshape(permute(l_chain, l_permute), [l_size, bond_size]);
    end

    function [r_chain_site, bond_size] = get_r_chain(r_map, r_map_2)

        r_map_2 = create_map(r_map_2, obj.numopts);
        r_map = create_map(r_map, obj.numopts);

        r_tensors = cell(1, r_map.N - 1);

        last_index = 0;

        for i = r_map.N:-1:2
            index_set = [0, 0, 0, 0];

            index_set_mask = r_map.leg_list{i}(3:end) > 0;
            non_empty = find(index_set_mask);

            if size(non_empty, 2) == 1
                last_index = non_empty(1);
                index_set(last_index) = r_map.N - i + 1;
            else
                prev_index = mod(last_index + 2, 4);

                if non_empty(1) == prev_index
                    index_set(non_empty(1)) = r_map.N - (i);
                    index_set(non_empty(2)) = r_map.N - (i - 1);
                    last_index = non_empty(2);
                else
                    index_set(non_empty(2)) = r_map.N - (i);
                    index_set(non_empty(1)) = r_map.N - (i - 1);
                    last_index = non_empty(1);
                end

            end

            r_tensors{i - 1} = obj.PEPO_cell{index_set(1) + 1, index_set(2) + 1, index_set(3) + 1, index_set(4) + 1};
        end

        leg_list_copy = r_map_2.leg_list;
        r_chain = ncon(r_tensors, leg_list_copy);

        r_chain_size = size(r_chain);
        r_size = d^(2 * r_map_2.N);

        bond_size = numel(r_chain) / r_size;

        r_chain_new_size = [r_chain_size(1:2 * r_map_2.N), bond_size];

        r_chain = reshape(r_chain, r_chain_new_size);

        r_permute = [2 * r_map_2.N + 1; site_ordering_permute(r_map_2.N)];
        r_chain_site = reshape(permute(r_chain, r_permute), [bond_size, r_size]);
    end

    function [x1, x2] = renumber(x, central, renumber)
        mask = x ~= 0;

        mask_central = x == central;

        x_size = size(x, 1);
        y = reshape(x(mask), [], 1);

        [~, y] = sort(y);
        y = reshape(y, x_size, []);

        x1 = x;
        x1(mask) = y;
        x2 = x1;
        x2(mask_central) = 0;

        if renumber == 1
            x2(~mask_central) = x2(~mask_central) - 1;
        end

    end

    if obj.testing == 1
        tensorarr = cell(1, 4);

        for ii = 1:4
            tensorarr{ii} = [1];
        end

    end

    for leg = 1:4
        inv_map = inv_maps{leg};

        if size(inv_map, 1) ~= 0

            if leg < 3
                central = max(inv_map);

                [i_map, i_map_2] = renumber(inv_map, central, 0);
                [l_chain_site, bond_size] = get_l_chain(i_map, i_map_2);

                if obj.testing == 1
                    tensorarr{leg} = l_chain_site;
                end

                if leg == 1
                    Tensor_site = reshape(Tensor_site, vector_sizes(1), []);
                    Tensor_site2 = l_chain_site \ Tensor_site;
                    vector_sizes(1) = bond_size;
                    Tensor_site = reshape(Tensor_site2, vector_sizes);
                else
                    Tensor_site = reshape(permute(Tensor_site, [2, 1, 3, 4, 5]), vector_sizes(2), []);
                    Tensor_site = l_chain_site \ Tensor_site;
                    vector_sizes(2) = bond_size;
                    Tensor_site = permute(reshape(Tensor_site, [vector_sizes(2), vector_sizes(1), vector_sizes(3), vector_sizes(4), vector_sizes(5)]), [2, 1, 3, 4, 5]);
                end

            else
                central = min(inv_map(inv_map ~= 0));

                [i_map, i_map_2] = renumber(inv_map, central, 1);
                [r_chain_site, bond_size] = get_r_chain(i_map, i_map_2);

                if obj.testing == 1
                    tensorarr{leg} = r_chain_site;
                end

                if leg == 4
                    Tensor_site = reshape(Tensor_site, [], vector_sizes(5));
                    Tensor_site = Tensor_site / r_chain_site;
                    vector_sizes(5) = bond_size;
                    Tensor_site = reshape(Tensor_site, vector_sizes);
                else
                    Tensor_site = reshape(permute(Tensor_site, [1, 2, 3, 5, 4]), [], vector_sizes(4));
                    Tensor_site2 = Tensor_site / r_chain_site;
                    vector_sizes(4) = bond_size;
                    Tensor_site = permute(reshape(Tensor_site2, [vector_sizes(1), vector_sizes(2), vector_sizes(3), vector_sizes(5), vector_sizes(4)]), [1, 2, 3, 5, 4]);
                end

            end

        end

    end

    if obj.testing
        Z = ncon({tensorarr{1}, tensorarr{2}, Tensor_site, tensorarr{3}, tensorarr{4}}, {[-1, 1], [-2, 2], [1, 2, -3, 3, 4], [3, -4], [4, -5]});

        ZZ = Z - Tensor_site_cpy;
    end

    if perm
        part = reshape(permute(Tensor_site, [3, 1, 2, 4, 5]), [d, d, vector_sizes(1), vector_sizes(2), vector_sizes(4), vector_sizes(5)]);
    else
        part = Tensor_site;
    end

end
