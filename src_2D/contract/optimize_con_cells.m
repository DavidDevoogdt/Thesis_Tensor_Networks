function [con_cells_cell2, targets] = optimize_con_cells(obj, maps, con_cells_cell, patterns, targets, ln_prefactor)
   


    if nargin < 6
        ln_prefactor = 0;
    end

    num_sub_probs = size(maps, 2);

    con_cells_cell2 = cell(1, num_sub_probs);

    for sub_prob = 1:num_sub_probs
        con_cells = con_cells_cell{sub_prob};
        con_cells2 = cell(1, 1);
        con_cells2_counter = 1;

        target = targets{sub_prob};
        map = maps{sub_prob};

        num_x = size(patterns, 2);

        for con_cell_index = 1:numel(con_cells)
            legs = con_cells{con_cell_index}{1};

            has_matched_pattern = 0;

            for pat_num = 1:num_x
                for ii = 1:size(legs, 2)
                    if same_pattern(legs{ii}, patterns{pat_num})
                        has_matched_pattern = 1;
                        break;
                    end
                end

                if has_matched_pattern
                    break;
                end
            end

            if has_matched_pattern %keep in new list
                con_cells2{con_cells2_counter} = con_cells{con_cell_index};
                con_cells2_counter = con_cells2_counter + 1;
            else %remove and change target
                temp_list_1 = fetch_PEPO_cells(obj, map, legs, ln_prefactor);

                A1 = ncon(temp_list_1, map.leg_list);
                perm_vect = [site_ordering_permute(map.N2); (2 * map.N2 + 1:size(size(A1), 2))'];

                target = target - reshape(permute(A1, perm_vect), size(target));
            end
        end

        con_cells_cell2{1, sub_prob} = con_cells2;
        targets{sub_prob} = target;
    end
end
