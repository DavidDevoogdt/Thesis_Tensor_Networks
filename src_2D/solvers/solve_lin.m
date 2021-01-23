function [x_cell, residual_target, rank_x,res_con] = solve_lin(obj, pattern, map, con_cells, target, lnprefact,loop_dim)
    
    if nargin < 6
        lnprefact = obj.nf;
    end

    %bring all parts without the PEPO cells to solve to the target
    [con_cells_cell2, target2] = optimize_con_cells(obj,{map}, {con_cells}, pattern, {target}, lnprefact);

    if numel(con_cells_cell2) ~= 1
        error("too many sub_problems")
    end

    if numel(con_cells_cell2{1}) ~= 1
        error("not linear")
    end

    residual_target = target2{1};
    cc = con_cells_cell2{1}{1};
    res_con =  con_cells_cell2{1};

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
        rem_map = remove_elem(num, rem_map);
    end
    
    [dims,dim_arr,bond_pairs] = removed_elems_dims(obj,num_pats,map,rem_map,pattern,nums);

   
    target_rot = permute_rhs(residual_target, nums); %put part in back

    function C = permute_rhs(B, nums)
        perm = 1:size(size(B), 2);
        perm(nums) = [];
        perm = [perm, nums'];

        C = permute(B, perm);
    end

    target_rot = reshape(target_rot, [], obj.dim^(2 * num_pats));
    dd = size(target_rot, 1);

    %cast problem to A_ij x_jk = b_ik and solve
    [A, ~] = contract_partial(obj,num, rem_map, {cc}, lnprefact);
    A_res = reshape(A, dd, []);
    x = lsqminnorm(A_res, target_rot, 10^(-12), 'nowarn');

    rank_x= rank(x);
    if rank_x == 0
        fprintf("rank 0 ")
    end

    % split into cells
    x = permute(reshape(x, [dim_arr, dimension_vector(obj.dim^2, num_pats)]), site_ordering_permute(num_pats, 1));
    x_cell = svd_x_cell(x,dims,bond_pairs,nums,loop_dim);

  
    
    
end
