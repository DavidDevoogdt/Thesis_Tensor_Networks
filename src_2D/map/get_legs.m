function leg_cell = get_legs(map, nums)
    %determine independent cycles for center tensor
    %used in solve_lin

    leg_cell = cell(0, 1);

    free = repmat([1, 2, 3, 4], numel(nums), 1);
    counter = 1;

    for row = 1:numel(nums)
        for col = 1:4

            %[row,col] = ind2sub( size(free),i);

            num = nums(row);

            involved = row;

            if free(row, col) == 0
                continue; %been here before
            end

            free(row, col) = 0;

            switch col
                case 1
                    n1 = map.h_bond_l_lookup{num};
                    if isempty(n1)
                        continue; %nothing to do
                    end
                    n1 = map.h_bonds{n1};
                case 2
                    n1 = map.v_bond_u_lookup{num};
                    if isempty(n1)
                        continue;
                    end
                    n1 = map.v_bonds{n1};
                case 3
                    n1 = map.h_bond_r_lookup{num};
                    if isempty(n1)
                        continue;
                    end
                    n1 = map.h_bonds{n1};
                case 4
                    n1 = map.v_bond_d_lookup{num};
                    if isempty(n1)
                        continue;
                    end
                    n1 = map.v_bonds{n1};
            end

            n1(n1 == num) = [];

            prev = num;
            curr = n1;

            if sum(nums == curr) ~= 0
                continue; %internal
            end

            c_nums = n1;

            while true
                neigh = next(prev, curr);
                switch numel(neigh)
                    case 0
                        break;
                    case 1
                        c_nums = [c_nums, neigh];
                        if neigh == c_nums
                            %fix involved, free, find bond
                            error('not implemented')
                            break;
                        end
                    otherwise
                        error('not implemented')
                end

                prev = curr;
                curr = neigh;
            end

            leg_cell{counter} = {{involved}, c_nums};

            counter = counter + 1;
        end
    end

    function neighours = next(prev, curr)
        neighours = [map.h_bonds{map.h_bond_l_lookup{curr}}, map.h_bonds{map.h_bond_r_lookup{curr}}, map.v_bonds{map.v_bond_u_lookup{curr}}, map.v_bonds{map.v_bond_d_lookup{curr}}];
        neighours(neighours == prev) = [];
        neighours(neighours == curr) = [];
    end
end
