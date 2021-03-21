function [F, G] = get_value_and_grad(obj, maps, con_cells_cell, root_patterns, extended_patterns, pattern_root, pattern_permutations, targets, x0, x_sizes, ln_prefactor)

    if nargin < 8
        ln_prefactor = 0;
    end

    num_sub_probs = size(maps, 2);

    F_cell = cell(num_sub_probs, 1);
    G_cell = cell(num_sub_probs, 1);

    F_total_size = 0;
    G_total_size = numel(x0);

    nargout_val = nargout;

    for sub_prob = 1:num_sub_probs
        con_cells = con_cells_cell{sub_prob};
        target = targets{sub_prob};
        map = maps{sub_prob};

        F_total_size = F_total_size + numel(target);

        num_x = size(x_sizes, 2);

        %x= reshape(x,[],1);

        x_cell = cell(num_x, 1);

        curr = 0;

        for i = 1:num_x
            num_elem = prod(x_sizes{i});
            x_cell{i} = reshape(x0(curr + 1:curr + num_elem), x_sizes{i});
            curr = curr + num_elem;
        end

        total_g_params = curr;

        F_sub_buffer = cell(1, size(con_cells, 2));
        G_sub_buffer = cell(1, size(con_cells, 2));

        %could be parforred
        for con_cell_index = 1:size(con_cells, 2)

            legs = con_cells{con_cell_index}{1};
            temp_list_1 = fetch_PEPO_cells(obj, map, legs, ln_prefactor, root_patterns, x_cell, extended_patterns, pattern_root, pattern_permutations);

            A1 = ncon(temp_list_1, map.leg_list);
            %F_sub = F_sub + reshape(  permute(A1,site_ordering_permute(map.N2)), size(target));

            perm_vect = [site_ordering_permute(map.N2); (2 * map.N2 + 1:size(size(A1), 2))'];

            F_sub_buffer{con_cell_index} = reshape(permute(A1, perm_vect), size(target));

            if nargout_val > 1 %calculate gradient

                G_sub_cell = cell(num_x, 1);

                for i = 1:num_x
                    num_elem = prod(x_sizes{i});
                    G_sub_cell{i} = zeros(numel(target), num_elem);
                end

                for pat_num = 1:num_x

                    x = x_cell{pat_num};

                    pat_num_mask = pattern_root == pat_num;
                    all_patterns = [{root_patterns{pat_num}}, extended_patterns{pat_num_mask}];
                    all_perm = [{1:4}, pattern_permutations{pat_num_mask}];

                    if size(root_patterns{pat_num}, 2) == 2 %boundary matrix

                        error("not verified for sub_num impelementation")

                        grad_total_size = [numel(target), numel(x)];
                        Grad_total = zeros(grad_total_size);

                        for ii = 1:size(legs, 2)

                            if same_pattern(legs{ii}, root_patterns{pat_num})

                                map2 = remove_elem(ii, map);

                                [Ai, ~] = contract_partial(obj, ii, map2, con_cells(con_cell_index), ln_prefactor, x_cell, root_patterns);

                                Grad_total = Grad_total + Ai;
                            end
                        end
                    else

                        size_x_red_o = size(x, [1, 2, 3, 4, 5, 6]);

                        d2 = obj.dim^2;
                        size_x_red_o(1) = d2;
                        size_x_red_o(2) = [];

                        grad_total_size = [size(target), size_x_red_o];
                        Grad_total = zeros(grad_total_size);

                        for sub_num = 1:numel(all_patterns)

                            sub_num_perm = [1, all_perm{sub_num} + 1];

                            size_x_red = size_x_red_o(sub_num_perm); %permute to current version

                            num = 0;

                            for ii = 1:size(legs, 2)
                                if same_pattern(legs{ii}, all_patterns{sub_num})

                                    index_before = num;
                                    index_after = map.N2 - num - 1;

                                    external_sizes = numel(target) / d2^(map.N2);

                                    map2 = remove_elem(ii, map);

                                    %[Ai, ~] = contract_partial(obj, ii, map2, con_cells(con_cell_index), ln_prefactor, x_cell, root_patterns);
                                    [Ai, ~] = contract_partial(obj, ii, map2, con_cells(con_cell_index), ln_prefactor, x_cell, root_patterns, extended_patterns, pattern_root, pattern_permutations);

                                    size_x = size_x_red(2:5);
                                    non_connected = map.leg_list{ii} < 0;
                                    non_connected = non_connected(3:end);
                                    size_x_internal = size_x;
                                    size_x_internal(non_connected) = 1;

                                    size_x_external = size_x;
                                    size_x_external(~non_connected) = 1;

                                    ai_size = size(Ai);
                                    external_size = ai_size(2:4);
                                    size_curr = external_sizes / prod(external_size);

                                    %permute everything to final form
                                    rvect_ext = [d2^index_before, d2, d2^index_after, external_size(1), size_curr, external_size(3), size_x_red];
                                    zz = numel(rvect_ext) - 4;
                                    rvect_ext = rvect_ext([1:zz, all_perm{sub_num} + zz]);

                                    Grad_total = reshape(Grad_total, rvect_ext);

                                    rvect = [d2^index_before, 1, d2^index_after, external_size, 1, size_x_internal];
                                    Ai_res = reshape(Ai, rvect);
                                    zz = numel(rvect) - 4;

                                    Ai_res = permute(Ai_res, [1:zz, all_perm{sub_num} + zz]);

                                    non_connected = non_connected(all_perm{sub_num});
                                    size_x_external = size_x_external(all_perm{sub_num});

                                    if isequal([0, 0, 0, 0], non_connected)

                                        for iii = 1:size_x_red(1)
                                            Grad_total(:, iii, :, :, 1, :, iii, :, :, :, :) = Grad_total(:, iii, :, :, 1, :, iii, :, :, :, :) ...
                                                + Ai_res(:, 1, :, :, 1, :, 1, :, :, :, :);
                                        end
                                    elseif isequal([0, 1, 1, 1], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [~, i2, i3, i4] = ind2sub(size_x_external, iiii);

                                                Grad_total(:, iii, :, :, iiii, :, iii, :, i2, i3, i4) = Grad_total(:, iii, :, :, iiii, :, iii, :, i2, i3, i4) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, :, 1, 1, 1);
                                            end
                                        end
                                    elseif isequal([1, 1, 0, 1], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [i1, i2, ~, i4] = ind2sub(size_x_external, iiii);
                                                Grad_total(:, iii, :, :, iiii, :, iii, i1, i2, :, i4) = Grad_total(:, iii, :, :, iiii, :, iii, i1, i2, :, i4) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, 1, 1, :, 1);
                                            end
                                        end
                                    elseif isequal([0, 1, 0, 1], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [~, i2, ~, i4] = ind2sub(size_x_external, iiii);

                                                Grad_total(:, iii, :, :, iiii, :, iii, :, i2, :, i4) = Grad_total(:, iii, :, :, iiii, :, iii, :, i2, :, i4) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, :, 1, :, 1);
                                            end
                                        end
                                    elseif isequal([0, 0, 1, 1], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [~, ~, i3, i4] = ind2sub(size_x_external, iiii);

                                                Grad_total(:, iii, :, :, iiii, :, iii, :, :, i3, i4) = Grad_total(:, iii, :, :, iiii, :, iii, :, :, i3, i4) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, :, :, 1, 1);
                                            end
                                        end
                                    elseif isequal([0, 1, 1, 0], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [~, i2, i3, ~] = ind2sub(size_x_external, iiii);

                                                Grad_total(:, iii, :, :, iiii, :, iii, :, i2, i3, :) = Grad_total(:, iii, :, :, iiii, :, iii, :, i2, i3, :) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, :, 1, 1, :);
                                            end
                                        end
                                    elseif isequal([1, 0, 0, 1], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [i1, ~, ~, i4] = ind2sub(size_x_external, iiii);

                                                Grad_total(:, iii, :, :, iiii, :, iii, i1, :, :, i4) = Grad_total(:, iii, :, :, iiii, :, iii, i1, :, :, i4) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, 1, :, :, 1);
                                            end
                                        end
                                    elseif isequal([1, 1, 0, 0], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [i1, i2, ~, ~] = ind2sub(size_x_external, iiii);

                                                Grad_total(:, iii, :, :, iiii, :, iii, i1, i2, :, :) = Grad_total(:, iii, :, :, iiii, :, iii, i1, i2, :, :) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, 1, 1, :, :);
                                            end
                                        end
                                    elseif isequal([0, 0, 0, 1], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [~, ~, ~, i4] = ind2sub(size_x_external, iiii);

                                                Grad_total(:, iii, :, :, iiii, :, iii, :, :, :, i4) = Grad_total(:, iii, :, :, iiii, :, iii, :, :, :, i4) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, :, :, :, 1);
                                            end
                                        end
                                    elseif isequal([0, 0, 1, 0], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [~, ~, i3, ~] = ind2sub(size_x_external, iiii);

                                                Grad_total(:, iii, :, :, iiii, :, iii, :, :, i3, :) = Grad_total(:, iii, :, :, iiii, :, iii, :, :, i3, :) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, :, :, 1, :);
                                            end
                                        end
                                    elseif isequal([0, 1, 0, 0], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [~, i2, ~, ~] = ind2sub(size_x_external, iiii);

                                                Grad_total(:, iii, :, :, iiii, :, iii, :, i2, :, :) = Grad_total(:, iii, :, :, iiii, :, iii, :, i2, :, :) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, :, 1, :, :);
                                            end
                                        end
                                    elseif isequal([1, 0, 0, 0], non_connected)
                                        for iii = 1:size_x_red(1)
                                            for iiii = 1:size_curr
                                                [i1, ~, ~, ~] = ind2sub(size_x_external, iiii);

                                                Grad_total(:, iii, :, :, iiii, :, iii, i1, :, :, :) = Grad_total(:, iii, :, :, iiii, :, iii, i1, :, :, :) ...
                                                    + Ai_res(:, 1, :, :, 1, :, 1, 1, :, :, :);
                                            end
                                        end
                                    else
                                        error("not implemented")

                                    end

                                end

                                Grad_total = reshape(Grad_total, grad_total_size);

                                if size(legs{ii}, 2) == 4 %not boundary matrix
                                    num = num + 1;
                                end

                            end
                        end
                    end

                    G_sub_cell{pat_num} = G_sub_cell{pat_num} +reshape(Grad_total, size(G_sub_cell{pat_num}));
                end

                G_sub_local = zeros([numel(target), total_g_params]);

                curr = 0;

                for i = 1:num_x
                    num_elem = prod(x_sizes{i});
                    G_sub_local(:, curr + 1:curr + num_elem) = G_sub_cell{i};
                    curr = curr + num_elem;

                end

                G_sub_buffer{con_cell_index} = G_sub_local;
            end
        end

        F_sub = -target;

        if nargout_val > 1
            G_sub = 0;
        end

        for con_cell_index = 1:size(con_cells, 2)
            F_sub = F_sub + F_sub_buffer{con_cell_index};

            if nargout_val > 1
                G_sub = G_sub + G_sub_buffer{con_cell_index};
            end
        end

        %put back toghether
        F_cell{sub_prob} = F_sub;

        if nargout_val > 1
            G_cell{sub_prob} = G_sub;
        end
    end

    F = zeros(F_total_size, 1);

    if nargout_val > 1
        G = zeros(F_total_size, G_total_size);
    end

    curr = 0;

    for sub_prob = 1:num_sub_probs
        num_elem = numel(F_cell{sub_prob});
        F(curr + 1:curr + num_elem) = reshape(F_cell{sub_prob}, [num_elem, 1]);

        if nargout > 1
            G(curr + 1:curr + num_elem, :) = G_cell{sub_prob};
        end

        curr = curr + num_elem;
    end
end
