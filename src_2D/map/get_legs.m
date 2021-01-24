function leg_cell  = get_legs(map,center)
    %determine independent cycles for center tensor
    %used in solve_lin

    leg_cell = cell(4,1);

    free = [1,2,3,4]

    counter = 1;
    while ~isempty(free)

        dir = free(1);
        free(1) = [];

        switch dir
            case 1
                n1 = map.h_bond_l_lookup{center};
                if isempty(n1)
                    leg_cell{counter} = {dir, [] };
                    counter = counter+1;
                    continue; 
                 end
                n1 = map.h_bonds{n1};
            case 2
                n1 = map.v_bond_u_lookup{center};
                if isempty(n1)
                    leg_cell{counter} = {dir, [] };
                    counter = counter+1;
                    continue; 
                 end
                n1 = map.v_bonds{n1};
            case 3 
                n1 = map.h_bond_r_lookup{center};
                if isempty(n1)
                    leg_cell{counter} = {dir, [] };
                    counter = counter+1;
                    continue; 
                 end
                n1 = map.h_bonds{n1};
            case 4
                n1 = map.v_bond_d_lookup{center};
                if isempty(n1)
                    leg_cell{counter} = {dir, [] };
                    counter = counter+1;
                    continue; 
                 end
                n1 = map.v_bonds{n1};
        end

        
        n1(n1==center) = [];

        prev = center;
        curr = n1;

        nums = n1;

        while true
            neigh = next(prev,curr);
            switch numel(neigh)
                case 0
                    break;
                case 1
                    nums = [nums,neigh];
                    if neigh==center
                        break;
                    end
                otherwise
                    error('not implemented')
            end

            prev = curr;
            curr = neigh;
        end

        leg_cell{counter} = {dir,nums};

        

        counter = counter+1;
    end

    function neighours =  next(prev,curr)
        neighours = [ map.h_bonds{map.h_bond_l_lookup{curr}},map.h_bonds{map.h_bond_r_lookup{curr}},map.v_bonds{map.v_bond_u_lookup{curr}},map.v_bonds{map.v_bond_d_lookup{curr}}];
        neighours( neighours==prev ) = [];
        neighours( neighours==curr) = [];
    end
end

