function [x_cell, residual_target, rank_x, res_con] = solve_lin(obj, pattern, map, con_cells, target, lnprefact, loop_dim, loop)

    if nargin < 6
        lnprefact = obj.nf;
    end

    if nargin < 8
        loop = 0;
    end

    %bring all parts without the PEPO cells to solve to the target
    [con_cells_cell2, target2] = optimize_con_cells(obj, {map}, {con_cells}, pattern, {target}, lnprefact);

    if numel(con_cells_cell2) ~= 1
        error("too many sub_problems")
    end

    if numel(con_cells_cell2{1}) ~= 1
        error("not linear")
    end

    residual_target = target2{1};
    cc = con_cells_cell2{1}{1};
    res_con = con_cells_cell2{1};

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

    rem_map = map;

    for pat = 1:num_pats
        num = nums(pat);
        rem_map = remove_elem(num, rem_map);
    end

    [dims, dim_arr, bond_pairs, ext_dims] = removed_elems_dims(obj, num_pats, map, rem_map, pattern, nums);

    method = nargin < 6;

    inv_eps = obj.inv_eps;

    if loop == 0 %invert leg per leg

        %x_sol =  residual_target;

        tensors = fetch_PEPO_cells(obj, map, cc{1}, lnprefact);

        %legs2 = get_legs(map, nums);

        legs = get_legs_2(rem_map, nums);
        %

        if ~isempty(legs)

            ext_ord2 = [];

            perm_vect = zeros(map.N - num_pats, 1);
            perm_dims = zeros(size(legs, 1), 1);
            perm_c = 1;
            for ii = 1:numel(legs)

                ext_ord2 = [ext_ord2, legs{ii}{1}];

                leg = sort(legs{ii}{2});
                nl = size(leg, 2);
                perm_vect(perm_c:perm_c + nl - 1) = leg;
                perm_c = perm_c + nl;
                perm_dims(ii) = obj.dim^(2 * nl);
            end

            [~, order] = sort(-ext_ord2);
            eord = cell(numel(legs), 1);

            perm_dims = [perm_dims, obj.dim^(2 * num_pats)];
            perm_vect = [perm_vect; nums];

            %x_sol = reshape( permute( x_sol,  [perm_vect',  nums(1) ] ),  [perm_dims', obj.dim^2]);
            %x_sol_dims = size(x_sol);
            %perm_basis = 1:5;

            A_list = cell(numel(legs), 1);

            for ii = 1:numel(legs)
                leg = legs{ii}{2};

                mask = zeros(size(map.num_map));
                for iii = 1:numel(leg)
                    mask = mask + (map.num_map == leg(iii));
                end
                mask = mask ~= 0;

                new_map = map.num_map .* mask;

                %make newe map with only branch

                sz = size(new_map);
                new_map = reshape(new_map, [], 1);
                [~, ord] = sort(new_map(new_map > 0));

                new_map(new_map > 0) = iorder(ord);

                new_map = reshape(new_map, sz);

                new_map = create_map(new_map, obj.numopts);

                A = ncon(tensors(sort(leg)), new_map.leg_list);

                size_A = size(A);
                ext = size_A(2 * numel(leg) + 1:end);

                %ext_perm = find(ext ~= 1);
                d2 = [2 * numel(leg) + 1:size(size_A, 2)];
                %patch = d2(ext_perm);
                %d2(ext_perm) = patch(iorder(external_order_orig)  ) ;

                pvect = [site_ordering_permute(numel(leg))', d2];

                A = reshape(permute(A, pvect), ...
                    [obj.dim^(2 * numel(leg)), prod(ext)]);

                A_list{ii} = A;
            end

            target_rot = permute(residual_target, perm_vect);
            target_rot = reshape(target_rot, perm_dims);

            x = lin_solver_core(A_list, target_rot, inv_eps);

            %fix order

            ext_ord2 = [];
            for ii = 1:numel(legs)

                ext_ord2 = [ext_ord2, legs{ii}{1}];
            end

            [~, order] = sort(-ext_ord2);

            x = reshape(x, [ext_dims(iorder(order)), obj.dim^(2 * num_pats)]); %split in dims per leg

            x = permute(x, [order, numel(order) + 1]); %put legs in order according to sites

            x = reshape(x, [dim_arr, dimension_vector(obj.dim^2, num_pats)]); %put back per site

            x = permute(x, site_ordering_permute(num_pats, 1)); %put physical dims next to site

        else
            x = residual_target;
            x = reshape(x, [dim_arr, dimension_vector(obj.dim^2, num_pats)]);
            x = permute(x, site_ordering_permute(num_pats, 1));
        end
    else % do it all at once, slow for large cells
        %remove target from map to the back and adapt order of target

        target_rot = permute_rhs(residual_target, nums); %put part in back

        target_rot = reshape(target_rot, [], obj.dim^(2 * num_pats));
        dd = size(target_rot, 1);

        %cast problem to A_ij x_jk = b_ik and solve
        [A, ~] = contract_partial(obj, num, rem_map, {cc}, lnprefact);
        A_res = reshape(A, dd, []);
        %works always,but potentially slow

        dA = decomposition(A_res, 'qr', 'RankTolerance', inv_eps, 'CheckCondition', false);
        x = dA \ target_rot;

        x = reshape(x, [dim_arr, dimension_vector(obj.dim^2, num_pats)]);
        x = permute(x, site_ordering_permute(num_pats, 1));

        %x = lsqminnorm(A_res,target_rot,inv_eps);
    end

    rank_x = NaN;

    x_cell = svd_x_cell(x, dims, bond_pairs, nums, loop_dim);

end

function C = permute_rhs(B, nums)
    perm = 1:size(size(B), 2);
    perm(nums) = [];
    perm = [perm, nums'];

    C = permute(B, perm);
end
