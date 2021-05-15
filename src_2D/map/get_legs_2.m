function l = get_legs_2(map, nums)

    m_copy = reshape(map.num_map, [], 1);

    l = cell(0, 1);

    while (true)
        [m, idx] = max(m_copy);

        if m ~= 0

            list = [];
            ext_legs = {};

            points = {idx};

            while numel(points) ~= 0
                p = points{1};
                points(1) = [];

                curr = m_copy(p);

                if curr~=0
                
                    m_copy(p) = 0;

                    list = [list, curr];

                    el = map.leg_list{curr};
                    ext_legs{end + 1} = el(el <- map.external_orig);

                    o = n(curr);

                    a = find(o);
                    for j = 1:numel(a)
                        np = o(a(j));

                        idx = find(m_copy == np);
                        if ~isempty(idx)
                            points{end + 1} = idx;
                        end

                    end
                end

            end

            l{end + 1} = {ext_legs, list};

        else
            break;
        end

    end

    function neighours = next(prev, curr)
        neighours = [map.h_bonds{map.h_bond_l_lookup{curr}}, map.h_bonds{map.h_bond_r_lookup{curr}}, map.v_bonds{map.v_bond_u_lookup{curr}}, map.v_bonds{map.v_bond_d_lookup{curr}}];
        neighours(neighours == prev) = [];
        neighours(neighours == curr) = [];

    end

    function o = n(curr)
        m2 = {map.h_bond_l_lookup{curr}, map.v_bond_u_lookup{curr}, map.h_bond_r_lookup{curr}, map.v_bond_d_lookup{curr}};
        mask = ~cellfun(@isempty, m2);

        n = next(-1, curr);

        o = zeros(4, 1);
        o(mask) = n;

        for i = 1:numel(nums)
            o(o == nums(i)) = [];
        end

    end

    %reverse order of elements

    for ii = 1:numel(l)
        [a, idx] = sort(l{ii}{2});
        b = l{ii}{1};

        l{ii}{2} = a;
        b = b(idx);

        c = [];
        for iii = 1:numel(b)
            c = [c, b{iii}];
        end

        l{ii}{1} = c;

    end

end
