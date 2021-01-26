function map2 = remove_elem(num, map)
    %removes element from map and puts indices in the back
    %original contracted indces are still contract when both
    %elements are removed

    map2 = map;

    if ~isfield(map, 'ext_h_bond_l_lookup')
        map2.ext_h_bond_l_lookup = cell(map.N, 1);
        map2.ext_h_bond_r_lookup = cell(map.N, 1);
        map2.ext_v_bond_u_lookup = cell(map.N, 1);
        map2.ext_v_bond_d_lookup = cell(map.N, 1);
    end

    if ~isfield(map, 'external_orig')
        external_orig = map.external_legs;
        map2.external_orig = external_orig;
    else
        external_orig = map.external_orig;
    end

    if ~isfield(map, 'final_order')
        max_occupied = map.external_legs;
    else
        max_occupied = max(-map.final_order);
    end

    l = map.h_bond_l_lookup{num};
    r = map.h_bond_r_lookup{num};
    u = map.v_bond_u_lookup{num};
    d = map.v_bond_d_lookup{num};

    con_list_cpy = map.leg_list;

    ii = 0;

    %h_bonds_missing = [];
    %v_bonds_missing = [];

    if ~isempty(l)
        ii = ii + 1;
        pair = map.h_bonds{l};

        map2.h_bond_l_lookup{num} = [];
        map2.num_h_bonds = map2.num_h_bonds - 1;
        %h_bonds_missing = [h_bonds_missing, l];

        other = pair(1);

        map2.h_bond_r_lookup{other} = [];

        if map.is_x_border(other)
            con_list_cpy{other}(2) = -(max_occupied +ii);
        else
            con_list_cpy{other}(2 + 3) = -(max_occupied +ii);
        end
    end

    if ~isempty(u)
        ii = ii + 1;
        pair = map.v_bonds{u};

        map2.v_bond_u_lookup{num} = [];
        map2.num_v_bonds = map2.num_v_bonds - 1;
        %v_bonds_missing = [v_bonds_missing, u];

        other = pair(1);

        map2.v_bond_d_lookup{other} = [];

        if map.is_y_border(other)
            con_list_cpy{other}(2) = -(max_occupied +ii);
        else
            con_list_cpy{other}(2 + 4) = -(max_occupied +ii);
        end
    end

    if ~isempty(r)
        ii = ii + 1;
        pair = map.h_bonds{r};

        map2.h_bond_r_lookup{num} = [];
        map2.num_h_bonds = map2.num_h_bonds - 1;
        %h_bonds_missing = [h_bonds_missing, r];

        other = pair(2);

        map2.h_bond_l_lookup{other} = [];

        if map.is_x_border(other)
            con_list_cpy{other}(1) = -(max_occupied +ii);
        else
            con_list_cpy{other}(2 + 1) = -(max_occupied +ii);
        end
    end

    if ~isempty(d)
        ii = ii + 1;
        pair = map.v_bonds{d};

        map2.v_bond_u_lookup{num} = [];
        map2.num_v_bonds = map2.num_v_bonds - 1;
        %h_bonds_missing = [h_bonds_missing, d];

        other = pair(2);

        map2.v_bond_d_lookup{other} = [];

        if map.is_y_border(other)
            con_list_cpy{other}(1) = -(max_occupied +ii);
        else
            con_list_cpy{other}(2 + 2) = -(max_occupied +ii);
        end

    end

    leg_list_num = map.leg_list{num}; %check whether it has a bond with already removed element

    for i = 1:4

        if leg_list_num(2 + i) < -external_orig

            switch i
                case 1
                    map2.ext_h_bond_l_lookup{num} = true;
                case 2
                    map2.ext_v_bond_u_lookup{num} = true;
                case 3
                    map2.ext_h_bond_r_lookup{num} = true;
                case 4
                    map2.ext_v_bond_d_lookup{num} = true;
            end
        end
    end

    num_legs_cpy = con_list_cpy{num};

    if ~isfield(map, 'leg_list_mask')
        map2.leg_list_mask = ones(map.N, 1) == 1;
        map2.leg_list_mask(num) = 0;
    else
        map2.leg_list_mask = map.leg_list_mask;
        map2.leg_list_mask(num) = 0;
    end

    %x0_list = con_list_cpy{ num};
    %con_list_cpy( num) = [];

    if ~isfield(map, 'final_external_missing')
        final_external_missing = [];
    else
        final_external_missing = map.final_external_missing;
    end

    if ~isfield(map, 'final_internal_missing')
        final_internal_missing = [];
    else
        final_internal_missing = map.final_internal_missing;
    end

    final_order = -1:-1:-(max_occupied + ii);
    seq = 1:(map.internal_legs + numel(final_internal_missing)); %contraction sequence

    final_external_missing = sort([final_external_missing, - num_legs_cpy(num_legs_cpy < 0)], 'descend');
    final_internal_missing = sort([final_internal_missing, num_legs_cpy(num_legs_cpy > 0)], 'descend');

    for l = 1:size(final_external_missing, 2)
        final_order(final_external_missing(l)) = [];
    end

    for l = 1:size(final_internal_missing, 2)
        seq(final_internal_missing(l)) = [];
    end

    %repackage in a new map

    map2.internal_legs = map2.internal_legs - ii;
    map2.external_legs = map2.external_legs - 4 + ii;

    map2.h_bond_l_lookup{num} = [];
    map2.h_bond_r_lookup{num} = [];
    map2.v_bond_u_lookup{num} = [];
    map2.v_bond_d_lookup{num} = [];

    map2.final_external_missing = final_external_missing;
    map2.final_internal_missing = final_internal_missing;

    map2.final_order = final_order;
    map2.seq = seq;

    %map2.h_bonds(h_bonds_missing) = [0,0];
    %map2.v_bonds(v_bonds_missing) = [];

    map2.leg_list = con_list_cpy;

    if ~isfield(map, "ii")
        map2.ii = ii;
    else
        map2.ii = map2.ii + ii;
    end
end
