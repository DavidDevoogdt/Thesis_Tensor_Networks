function [dims, dim_arr, bond_pairs,ext_dims] = removed_elems_dims(obj, num_pats, map, rem_map, pattern, nums)

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
    
    ext_dims = [];
    
    
    for i=1:num_pats
        
        
        
        d2 = dims{i}(3:end);
        
        mask =  d2~=-1;
        mask2 = d2~=1;
        mask3 =mask & mask2;
        
        ext_dims = [ext_dims,  d2(mask3) ];
        
    end
    
end
