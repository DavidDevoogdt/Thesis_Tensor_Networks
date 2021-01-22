function fixed_bonds = get_fixed_bonds(map, bonds)
    %bonds is nx1 cell with 2x1 cells with {[n1,n2], virt_level}

    num_fixed = size(bonds, 2);

    fixed_bonds = zeros(map.internal_legs, 1) - 1;

    for i = 1:num_fixed
        bond = bonds{i};
        n = bond{1};
        n1 = n(1); n2 = n(2);
        virt_level = bond{2};

        [h_bond, ~] = intersect([map.h_bond_r_lookup{n1}, map.h_bond_l_lookup{n1}], [map.h_bond_r_lookup{n2}, map.h_bond_l_lookup{n2}]);
        [v_bond, ~] = intersect([map.v_bond_u_lookup{n1}, map.v_bond_d_lookup{n1}], [map.v_bond_u_lookup{n2}, map.v_bond_d_lookup{n2}]);

        for ii = 1:size(h_bond, 1)
            fixed_bonds(h_bond(ii)) = virt_level;
        end

        for ii = 1:size(v_bond, 1)
            fixed_bonds(v_bond(ii) + map.num_h_bonds) = virt_level;
        end

        pot_bond_v{1} = map.v_bond_u_lookup{n1};
        pot_bond_v{2} = map.v_bond_u_lookup{n2};
    end
end
