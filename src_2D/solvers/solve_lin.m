function [x_cell, residual_target, rank_x, mask] = solve_lin(obj, pattern, map, con_cells, target, lnprefact, disable_check)
    %x is a tensor with the solved part, dims are the individual
    %dimensions of the PEPO cells. The still connected bonds are
    %numbered with negative indices in dims

    if nargin < 6
        lnprefact = obj.nf;
    end

    if nargin < 7
        disable_check = false;
    end

    %bring all parts without the PEPO cells to solve to the target
    [con_cells_cell2, target2] = optimize_con_cells(obj,{map}, {con_cells}, pattern, {target}, lnprefact);

    if numel(con_cells_cell2) ~= 1
        error("to many sub_problems")
    end

    if numel(con_cells_cell2{1}) ~= 1
        error("not linear")
    end

    residual_target = target2{1};

    cc = con_cells_cell2{1}{1};

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

    %remove target from map to the back and adapt order of target
    rem_map = map;

    for pat = 1:num_pats
        num = nums(pat);
        rem_map = remove_elem(obj,num, rem_map);

    end

    target_rot = rotate_rhs(residual_target, nums); %put part in back

    function C = rotate_rhs(B, nums)
        perm = 1:size(size(B), 2);
        perm(nums) = [];
        perm = [perm, nums'];

        C = permute(B, perm);
    end

    %determine appropriate size of differnten PEPO cells and
    %connections between them
    dims = cell(1, num_pats);

    dim_arr = zeros(1, num_pats);

    for n1 = 1:num_pats
        p = pattern{n1};
        ll = map.leg_list{nums(n1)}(3:6);

        dd = 1;

        for j = 1:4

            if mod(j, 2) == 1
                d = obj.virtual_level_sizes_horiz(p(j) + 1);
            else
                d = obj.virtual_level_sizes_vert(p(j) + 1);
            end

            dims{n1} = [dims{n1}, d];
            dd = dd * d;
        end

        dims{n1} = [obj.dim, obj.dim, dims{n1}];

        dim_arr(n1) = dd;
    end

    bond_counter = 0;

    bond_pairs = {};

    for n1 = 1:num_pats%look for bonds between extracted x
        num = nums(n1);

        for j = 1:4

            if ll(j) > 0

                switch j
                    case 1

                        if ~isempty(rem_map.ext_h_bond_l_lookup{num})
                            bond_counter = bond_counter +1;

                            bond = map.h_bond_l_lookup{num};
                            pair = map.h_bonds{bond};
                            other = pair(pair ~= num);

                            n2 = find(nums == other);

                            dim1 = dim_arr(n1);
                            dim2 = dim_arr(n2);
                            jdim1 = dims{n1}(2 + 1);
                            jdim2 = dims{n2}(2 + 3);

                            dim_arr(n1) = dim1 / jdim1;
                            dim_arr(n2) = dim2 / jdim2;
                            dims{n1}(3) = -bond_counter;
                            dims{n2}(5) = -bond_counter;

                            bond_pairs{bond_counter} = pair;

                        end

                    case 2

                        if ~isempty(rem_map.ext_v_bond_u_lookup{num})
                            bond_counter = bond_counter +1;

                            bond = map.v_bond_u_lookup{num};
                            pair = map.v_bonds{bond};
                            other = pair(pair ~= num);

                            n2 = find(nums == other);

                            dim1 = dim_arr(n1);
                            dim2 = dim_arr(n2);
                            jdim1 = dims{n1}(2 + 2);
                            jdim2 = dims{n2}(2 + 4);

                            dim_arr(n1) = dim1 / jdim1;
                            dim_arr(n2) = dim2 / jdim2;
                            dims{n1}(2 + 2) = -bond_counter;
                            dims{n2}(2 + 4) = -bond_counter;

                            bond_pairs{bond_counter} = pair;

                        end

                    case 3

                        if ~isempty(rem_map.ext_h_bond_r_lookup{num})
                            bond_counter = bond_counter +1;

                            bond = map.h_bond_r_lookup{num};
                            pair = map.h_bonds{bond};
                            other = pair(pair ~= num);

                            n2 = find(nums == other);

                            dim1 = dim_arr(n1);
                            dim2 = dim_arr(n2);
                            jdim1 = dims{n1}(2 + 3);
                            jdim2 = dims{n2}(2 + 1);

                            dim_arr(n1) = dim1 / jdim1;
                            dim_arr(n2) = dim2 / jdim2;
                            dims{n1}(2 + 3) = -bond_counter;
                            dims{n2}(2 + 1) = -bond_counter;

                            bond_pairs{bond_counter} = pair;

                        end

                    case 4

                        if ~isempty(rem_map.ext_v_bond_d_lookup{num})
                            bond_counter = bond_counter +1;

                            bond = map.v_bond_d_lookup{num};
                            pair = map.v_bonds{bond};
                            other = pair(pair ~= num);

                            n2 = find(nums == other);

                            dim1 = dim_arr(n1);
                            dim2 = dim_arr(n2);
                            jdim1 = dims{n1}(2 + 4);
                            jdim2 = dims{n2}(2 + 2);

                            dim_arr(n1) = dim1 / jdim1;
                            dim_arr(n2) = dim2 / jdim2;
                            dims{n1}(2 + 4) = -bond_counter;
                            dims{n2}(2 + 2) = -bond_counter;

                            bond_pairs{bond_counter} = pair;

                        end

                end

            end

        end

    end

    perm_vect = site_ordering_permute(num_pats, 1);

    target_rot = reshape(target_rot, [], obj.dim^(2 * num_pats));
    dd = size(target_rot, 1);

    %cast to A*x=b
    [A, ~] = contract_partial(obj,num, rem_map, {cc}, lnprefact);
    A_res = reshape(A, dd, []);

    %time intensive step

    ext_dim = size(target_rot, 2);

    %t_rot_norms = svds(target_rot,1);

    t_rot_norms = max(target_rot .* conj(target_rot), [], 1).^0.5;

    %             tr2 = reshape( target_rot,dimension_vector(obj.dim,2*map.N)  );
    %             tr2 = ipermute( tr2, site_ordering_permute(map.N) );
    %             tr2 = reshape(tr2, [obj.dim^map.N , obj.dim^map.N]);
    %
    %             tr_svds_norm = svds(tr2,1);

    A_res_norm = svds(A_res, 1);

    %max_x_norm = t_rot_norms/A_res_norm;

    n = 13; %

    n_min = 10;

    convergence_factor = exp((lnprefact - obj.nf) * num_pats);

    criterium = ones(1, ext_dim); %*convergence_factor;

    mask = ~zeros(1, ext_dim);

    x = lsqminnorm(A_res, target_rot, 10^(-12), 'nowarn');

    %           x = zeros( size(A_res,2), size(target_rot,2) );
    %             while n>n_min
    %                 %y = lsqminnorm(A_res,target_rot(:,mask),10^(-n),'nowarn');
    %                 y = lsqminnorm(A_res,target_rot,10^(-n),'nowarn');
    %                 %ny= sum(y.*conj(y), 1).^0.5;
    %
    %                 ny = max(abs(y));
    %
    %                 %y2 = reshape( y,dimension_vector(obj.dim,2*map.N)  );
    %                 %y2 = ipermute( y2, site_ordering_permute(map.N) );
    %                 %y2 = reshape(y2, [obj.dim^map.N , obj.dim^map.N]);
    %                 %ny = svds(y2,1 );
    %
    %
    %                 alpha = ny<criterium(mask);% criterium(mask);
    %
    %
    %                 %new_mask = mask;
    %                 %new_mask(new_mask==1)= alpha;
    %
    %                 %x(:,new_mask) = y(:,alpha);
    %
    %                 x=y;
    %
    %                 %mask(new_mask)=0;
    %
    %                 %if  sum(mask) == 0 || disable_check == true
    %                 if  sum(~alpha)== 0 || disable_check == true
    %
    %                    break;
    %                 else
    %                     fprintf(".")
    %                     n=n-1;
    %
    %                 end
    %
    %
    % %            end
    % %
    %             end
    %
    %
    %
    %             if ~(  sum(~alpha)== 0 || disable_check == true)
    %             %if  sum(mask) ~= 0 && disable_check == false
    %                     fprintf(",")
    %                     fprintf("%.4e,%.4e",max(ny),convergence_factor)
    %
    % %
    % %                     f_Arr = find(mask==1);
    % %
    % %                     for i = 1:size(f_Arr,2)
    % %                         ind = f_Arr(i);
    % %
    % %
    % %                         switch num_pats
    % %
    % %
    % %                         %fact =min(1/sqrt(ny(i)),sqrt(t_rot_norms(ind)/ny(i))*2/A_res_norm);
    % %
    % %
    % %                         x(:,ind) = y(:,i).*fact ;
    % %                     end
    % %
    % %                     %do step in direction of x st norm is still small
    % %
    %              end

    rank_x = rank(x);

    if rank_x == 0
        fprintf("rank 0 ")
    end

    x = permute(reshape(x, [dim_arr, dimension_vector(obj.dim^2, num_pats)]), perm_vect);

    %now split across each bond with svd

    switch bond_counter
        case 0
            x_cell = {reshape(x, dims{1})};
        case 1
            pair = bond_pairs{1};

            num1 = pair(1);
            num2 = pair(2);

            i1 = find(nums == num1);
            i2 = find(nums == num2);

            dims1 = dims{i1};
            dims2 = dims{i2};

            mask1 = dims1 == -1;
            mask2 = dims2 == -1;

            d1 = prod(dims1(~mask1));
            d2 = prod(dims2(~mask2));

            assert(d1 == d2)

            n1 = find(mask1);
            n2 = find(mask2);

            dim1_alt = [prod(dims1(1:n1 - 1)), prod(dims1(n1 + 1:end))];
            dim2_alt = [prod(dims2(1:n2 - 1)), prod(dims2(n2 + 1:end))];

            [U, S, V] = svd(reshape(x, d1, d1));

            sqrt_S = diag(diag(S).^0.5);

            l = permute(reshape(U * sqrt_S, dim1_alt(1), dim1_alt(2), []), [1, 3, 2]);
            r = permute(reshape(sqrt_S * V', [], dim2_alt(1), dim2_alt(2)), [2, 1, 3]);

            dims1(mask1) = size(l, 2);
            dims2(mask2) = size(r, 2);

            x_cell{i1} = reshape(l, dims1);
            x_cell{i2} = reshape(r, dims2);
        otherwise
            error("not implemted")

    end

end
