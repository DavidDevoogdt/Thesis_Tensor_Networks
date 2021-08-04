function [x_cell, residual_target, res_con] = solve_lin(obj, pattern, map, con_cells, target, lnprefact, split_opts, inv_opts)
    %code to calculate fast pseudoinverses for a block. use solve_lin_and_assign

    if nargin < 6
        lnprefact = obj.nf;
    end

    if nargin < 8
        inv_opts.full_inverse = 0;
        inv_opts.sigma = 1e-14;
    end

    %bring all parts without the PEPO cells to solve to the target
    [con_cells_cell2, target2] = optimize_con_cells(obj, {map}, {con_cells}, pattern, {target}, lnprefact);

    if numel(con_cells_cell2) ~= 1
        error("too many maps")
    end

    if numel(con_cells_cell2{1}) ~= 1
        error("pattern occurs multiple times, try solve non lin")
    end

    residual_target = target2{1};
    cc = con_cells_cell2{1}{1};
    res_con = con_cells_cell2{1};

    if isempty(cc)
        error('no valid combinations with pattern, try other one')
    end

    %get locaton of patterns
    num_pats = size(pattern, 2);
    nums = zeros(num_pats, 1) -1;

    for pat = 1:num_pats
        for i = 1:map.N
            if same_pattern(cc{1}{i}, pattern{pat})
                nums(pat) = i;
                break;
            end
        end
    end

    if sum(nums == -1) ~= 0
        error("pattern not found")
    end

    %make map without the patterns
    rem_map = map;

    for pat = 1:num_pats
        num = nums(pat);
        rem_map = remove_elem(num, rem_map);
    end

    [dims, dim_arr, bond_pairs, ext_dims] = removed_elems_dims(obj, num_pats, map, rem_map, pattern, nums);

    if inv_opts.full_inverse == 0 %invert leg per leg, faster

        tensors = fetch_PEPO_cells(obj, map, cc{1}, lnprefact);

        parts = get_disconnected_parts(rem_map, nums);

        if ~isempty(parts) % 0--|--1--|--0 removes all blocks

            %determine and reoder indices according
            order_pattern_externsions = [];

            perm_vect = zeros(map.N - num_pats, 1);
            perm_dims = zeros(size(parts, 1), 1);
            perm_c = 1;
            for ii = 1:numel(parts)

                order_pattern_externsions = [order_pattern_externsions, parts{ii}{1}];

                leg = sort(parts{ii}{2});
                nl = size(leg, 2);
                perm_vect(perm_c:perm_c + nl - 1) = leg;
                perm_c = perm_c + nl;
                perm_dims(ii) = obj.dim^(2 * nl);
            end

            perm_dims = [perm_dims, obj.dim^(2 * num_pats)];
            perm_vect = [perm_vect; nums];

            %fetch the necesarry tensors per leg, and make a new map with only the given leg.
            A_list = cell(numel(parts), 1);

            for ii = 1:numel(parts)
                leg = parts{ii}{2};

                %make newe map with only branch
                mask = zeros(size(map.num_map));
                for iii = 1:numel(leg)
                    mask = mask + (map.num_map == leg(iii));
                end
                mask = mask ~= 0;

                new_map = map.num_map .* mask;
                sz = size(new_map);
                new_map = reshape(new_map, [], 1);
                [~, ord] = sort(new_map(new_map > 0));

                new_map(new_map > 0) = iorder(ord);
                new_map = reshape(new_map, sz);
                new_map = create_map(new_map, obj.numopts);

                %Contract leg
                A = ncon(tensors(sort(leg)), new_map.leg_list);

                %reorder: all physical site indices in fron, one index for connected to central tensor
                size_A = size(A);
                ext = size_A(2 * numel(leg) + 1:end); %not phsycial indices
                d2 = [2 * numel(leg) + 1:size(size_A, 2)];
                pvect = [site_ordering_permute(numel(leg))', d2];
                A = reshape(permute(A, pvect), ...
                    [obj.dim^(2 * numel(leg)), prod(ext)]);

                %
                A_list{ii} = A;
            end

            %transform target according to legs, and reshape
            target_rot = permute(residual_target, perm_vect);
            target_rot = reshape(target_rot, perm_dims);

            %perfor inversion
            x = lin_solver_core(A_list, target_rot, inv_opts.inv_eps);

            %reorder external indices
            [~, order_idx] = sort(-order_pattern_externsions);
            x = reshape(x, [ext_dims(iorder(order_idx)), obj.dim^(2 * num_pats)]); %split according to dims of legs, and physical dims at end
            x = permute(x, [order_idx, numel(order_idx) + 1]); %reoder the exteral legs as numbered by rem_map
            x = reshape(x, [dim_arr, dimension_vector(obj.dim^2, num_pats)]); %split in dimension of external legs of a given site. One inverted leg can be split over multiple sites
            x = permute(x, site_ordering_permute(num_pats, 1)); %put physical dims next to site

        else % trivial case of 2x1 cell
            x = residual_target;
            x = reshape(x, [dim_arr, dimension_vector(obj.dim^2, num_pats)]);
            x = permute(x, site_ordering_permute(num_pats, 1));
        end
    else % do the inversion all at once, slow for large cells

        %remove target from map to the back and adapt order\ of target

        target_rot = permute_rhs(residual_target, nums); %put part in back
        target_rot = reshape(target_rot, [], obj.dim^(2 * num_pats));
        dd = size(target_rot, 1);

        %cast problem to A_ij x_jk = b_ik and solve
        [A, ~] = differentiate_x(obj, num, rem_map, cc, lnprefact);
        A_res = reshape(A, dd, []);
        %works always,but potentially slow

        dA = decomposition(A_res, 'qr', 'RankTolerance', inv_opts.inv_eps, 'CheckCondition', false);
        x = dA \ target_rot;

        x = reshape(x, [dim_arr, dimension_vector(obj.dim^2, num_pats)]);
        x = permute(x, site_ordering_permute(num_pats, 1));

    end

    x_cell = svd_x_cell(x, dims, bond_pairs, nums, split_opts);

end

function C = permute_rhs(B, nums)
    perm = 1:size(size(B), 2);
    perm(nums) = [];
    perm = [perm, nums'];

    C = permute(B, perm);
end
