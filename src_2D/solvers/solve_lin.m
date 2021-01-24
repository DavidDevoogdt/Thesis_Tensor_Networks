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


    if num_pats == 1 %invert leg per leg

        x_sol =  residual_target;

        tensors = fetch_PEPO_cells(obj, map, cc{1}, lnprefact);
        
        
        legs = get_legs(map,nums(1));

        perm_vect = zeros(map.N - num_pats,1);
        perm_dims = zeros(size(legs,1),1);
        perm_c = 1;
        for ii = 1:size(legs,1)
            leg = legs{ii}{2};
            nl = size(leg,2); 
            perm_vect(perm_c:perm_c+nl-1) = leg;
            perm_c = perm_c + nl;
            perm_dims(ii) = obj.dim^(2*nl);
        end

        x_sol = reshape( permute( x_sol,  [perm_vect',  nums(1) ] ),  [perm_dims', obj.dim^2]); 
        x_sol_dims = size(x_sol);
        perm_basis = 1:5;

        for ii = 1:size(legs,1)
            leg = legs{ii}{2};

            if ~ isempty(leg)
                mask = zeros( size(map.num_map));
                for iii = 1:numel(leg)
                    mask = mask + ( map.num_map == leg(iii)) ;
                end
                mask = mask~=0;

                new_map = map.num_map.*mask;
                [~, new_map(mask)  ] =  sort( new_map(mask) ) ;


                new_map = create_map( new_map, obj.numopts);

                A = ncon( tensors(  sort(leg) ), new_map.leg_list  );

                size_A = size(A);
                ext = size_A( 2*numel(leg)+1:end);

                A =  reshape( permute(    A, [site_ordering_permute( numel(leg)  )', 2*numel(leg)+1: size(size_A,2)    ]),...
                                    [  obj.dim^(2*numel(leg) ), prod( ext)  ]  );


                dA = decomposition(A,'cod','RankTolerance',1e-12,'CheckCondition',false);

                perm_basis_2 = perm_basis;
                perm_basis_2(ii) = [];

                %if ii < 3

                    perm_basis_2 = [ii,perm_basis_2];

                    x_sol = permute(x_sol , perm_basis_2);
                    new_xsol_dims = size(x_sol);
                    x_sol = reshape( x_sol, x_sol_dims(ii),[]);

                    x_sol = dA\x_sol;
%                 else
% 
%                     perm_basis_2 = [perm_basis_2,ii];
% 
%                     x_sol = permute(x_sol , perm_basis_2);
%                     new_xsol_dims = size(x_sol);
%                     x_sol = reshape( x_sol, [],x_sol_dims(ii));
% 
%                     x_sol = x_sol/dA;
%                 end

                x_sol = ipermute( reshape( x_sol, new_xsol_dims  ),perm_basis_2);
            end

        end


        x = reshape( permute( x_sol  , [5,1,2,3,4]), [obj.dim,obj.dim,  size(x_sol, 1:4)  ])  ;
        x_cell = {x};

        rank_x = NaN;

    else% do it all at once
            %remove target from map to the back and adapt order of target
        rem_map = map;

        for pat = 1:num_pats
            num = nums(pat);
            rem_map = remove_elem(num, rem_map);
        end
        
        [dims,dim_arr,bond_pairs] = removed_elems_dims(obj,num_pats,map,rem_map,pattern,nums);


        target_rot = permute_rhs(residual_target, nums); %put part in back

       
  
        target_rot = reshape(target_rot, [], obj.dim^(2 * num_pats));
        dd = size(target_rot, 1);

        %cast problem to A_ij x_jk = b_ik and solve
        [A, ~] = contract_partial(obj,num, rem_map, {cc}, lnprefact);
        A_res = reshape(A, dd, []);
        %works always,but potentially slow
        dA = decomposition(A_res,'cod','RankTolerance',1e-12,'CheckCondition',false);
        x = dA\target_rot;
        
        rank_x= rank(x);
        if rank_x == 0
            fprintf("rank 0 ")
        end

        % split into cells
        x = permute(reshape(x, [dim_arr, dimension_vector(obj.dim^2, num_pats)]), site_ordering_permute(num_pats, 1));
        x_cell = svd_x_cell(x,dims,bond_pairs,nums,loop_dim);
    end



end

function C = permute_rhs(B, nums)
    perm = 1:size(size(B), 2);
    perm(nums) = [];
    perm = [perm, nums'];

    C =  permute(B, perm);
end

