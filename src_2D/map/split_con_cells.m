function [con_cells_1, con_cells_2] = split_con_cells(map, con_cells)

    con_cells_1 = {};
    con_cells_2 = {};

    for i = 1:size(con_cells, 2)
        con_cell = con_cells{i}{1};
        n_cells = size(con_cell, 2);
        neighbour = zeros(n_cells, 1);

        for j = 1:n_cells
            neighbour(j) = sum(con_cell{j} ~= 0);
        end

        end_points = find(neighbour == 1);

        covered = neighbour;

        for j = 1:size(end_points, 1)
            curr_point = end_points(j);
            n = 1;

            bool = 1;

            while bool
                cell = con_cell{curr_point};

                loc = find(cell == n);

                if isempty(loc)

                    bool = 0;

                else
                    covered(curr_point) = 0;

                    switch loc
                        case 1
                            bond = map.h_bond_l_lookup{curr_point};
                            other = map.h_bonds{bond};
                        case 2
                            bond = map.v_bond_u_lookup{curr_point};
                            other = map.v_bonds{bond};
                        case 3
                            bond = map.h_bond_r_lookup{curr_point};
                            other = map.h_bonds{bond};
                        case 4
                            bond = map.v_bond_d_lookup{curr_point};
                            other = map.v_bonds{bond};
                    end
                    curr_point = other(other ~= curr_point);
                end

                n = n + 1;
            end
        end

        good_cell = sum(covered ~= 0) <= 1;

        if good_cell
            con_cells_1{end + 1} = con_cells{i};
        else
            con_cells_2{end + 1} = con_cells{i};
        end
    end
end
