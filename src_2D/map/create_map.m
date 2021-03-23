function [map, boundary_map] = create_map(pos_map, opts, internal)
    %h_cyclic -> right column is boundary matrix mx
    %V_cyclic -> under row is boundary matrix my
    %matrices on the boundary should recieve the higher number than
    %the others

    if nargin < 2
        opts.numbered = 0;
        opts.v_cyclic = 0;
        opts.h_cyclic = 0;
        opts.boundary_matrix = 0;
    end

    if nargin < 3
        internal = 0; %recursive call
    end

    if ~isfield(opts, 'numbered')
        opts.numbered = 0;
    end

    if ~isfield(opts, 'v_cyclic')
        opts.v_cyclic = 0;
    end

    if ~isfield(opts, 'h_cyclic')
        opts.h_cyclic = 0;
    end

    if ~isfield(opts, 'boundary_matrix')
        opts.boundary_matrix = 0;
    end

    %number the location of operators from up to down and left to
    %right, and create toghether with it a leg_list for ncon
    map.opts = opts;

    [m, n] = size(pos_map);
    map.m = m;
    map.n = n;

    map.pos_lookup = {};

    if opts.numbered == 0
        counter = 1;
        %boundary matrices have highest numbers
        for x = 1:n - 1
            for y = 1:m - 1

                if pos_map(y, x) == 1
                    pos_map(y, x) = counter;
                    map.pos_lookup{counter} = [y, x];
                    counter = counter +1;
                end
            end
        end

        for y = 1:m - 1
            if pos_map(y, n) == 1
                pos_map(y, n) = counter;
                map.pos_lookup{counter} = [y, n];
                counter = counter +1;
            end
        end

        for x = 1:n
            if pos_map(m, x) == 1
                pos_map(m, x) = counter;
                map.pos_lookup{counter} = [m, x];
                counter = counter +1;
            end
        end

        N = counter - 1;
    else %todo more checking
        N = 0;

        for x = 1:n
            for y = 1:m
                if pos_map(y, x) ~= 0
                    counter = pos_map(y, x);
                    %pos_map(y,x) = counter;
                    map.pos_lookup{counter} = [y, x];

                    if counter > N
                        N = counter;
                    end
                end
            end
        end
    end

    %N number of physical sites

    map.N = N;

    map.num_map = pos_map;

    map.h_bonds = {};
    map.v_bonds = {};

    %left right up down bonds
    map.h_bond_l_lookup = cell(N, 1);
    map.h_bond_r_lookup = cell(N, 1);
    map.v_bond_u_lookup = cell(N, 1);
    map.v_bond_d_lookup = cell(N, 1);

    map.is_x_border = zeros(map.N, 1);
    map.is_y_border = zeros(map.N, 1);

    leg_list = cell(1, N);
    leg_list(:) = {[0, 0, 0, 0, 0, 0]};

    %do the horizontal internal bonds

    internal_counter = 1;

    for num = 1:N
        coor = map.pos_lookup{num};
        x = coor(2); y = coor(1);

        if x == n
            if opts.h_cyclic == 1
                next_x = 1;
                map.is_x_border(n) = 1;
            else
                continue; %skip this round
            end
        else
            next_x = x + 1;
        end

        if pos_map(y, x) ~= 0 && pos_map(y, next_x) ~= 0
            n1 = pos_map(y, x);
            n2 = pos_map(y, next_x);

            leg_list{n1}(5) = internal_counter;
            leg_list{n2}(3) = internal_counter;

            map.h_bonds{internal_counter} = [n1, n2];

            %save bonds per site
            map.h_bond_r_lookup{n1} = [map.h_bond_r_lookup{n1}, internal_counter];
            map.h_bond_l_lookup{n2} = [map.h_bond_l_lookup{n2}, internal_counter];

            internal_counter = internal_counter + 1;

            %fprintf( "hor %d-%d\n",pos_map(y,x), pos_map(y,x+1));
        end
    end

    map.num_h_bonds = internal_counter - 1;

    %vertical internal bonds

    for num = 1:N
        coor = map.pos_lookup{num};
        x = coor(2); y = coor(1);

        if y == m
            if opts.v_cyclic == 1
                next_y = 1;
                map.is_y_border(n) = 1;
            else
                continue; %skip this round
            end
        else
            next_y = y + 1;
        end

        if pos_map(y, x) ~= 0 && pos_map(next_y, x) ~= 0
            n1 = pos_map(y, x);
            n2 = pos_map(next_y, x);

            leg_list{n1}(6) = internal_counter;
            leg_list{n2}(4) = internal_counter;

            v_number = internal_counter - map.num_h_bonds;
            map.v_bonds{v_number} = [n1, n2];

            %save bonds per site
            map.v_bond_d_lookup{n1} = [map.v_bond_d_lookup{n1}, v_number];
            map.v_bond_u_lookup{n2} = [map.v_bond_u_lookup{n2}, v_number];

            internal_counter = internal_counter + 1;

            %fprintf( "vert %d-%d\n",pos_map(y,x), pos_map(y+1,x));
        end
    end

    map.num_v_bonds = internal_counter - map.num_h_bonds - 1;
    map.internal_legs = internal_counter - 1;

    if internal == 0 && opts.boundary_matrix == 1

        total_borders = sum(map.is_x_border) + sum(map.is_y_border);
        map.N2 = map.N - total_borders;

        if sum(map.is_x_border(1:map.N2)) ~= 0
            error("boundary matrices should have largest numbers")
        end

        if sum(map.is_y_border(1:map.N2)) ~= 0
            error("boundary matrices should have largest numbers")
        end

    else %boundary row and column already deleted.
        map.N2 = map.N;
        map.is_x_border = map.is_x_border * 0;
        map.is_y_border = map.is_y_border * 0;
    end

    %number ij according to number
    external_counter = 1;

    for N1 = 1:map.N2

        leg_list{N1}(1) = -(N1);
        leg_list{N1}(2) = -(N1 + map.N2);

        external_counter = external_counter + 1;
    end

    for N1 = map.N2 + 1:map.N
        leg_list{N1} = leg_list{N1}(leg_list{N1} ~= 0);
    end

    external_counter = 2 * (map.N2) + 1;

    %number all other indices

    for i = 1:N
        if map.is_x_border(i) == 0 && map.is_y_border(i) == 0
            for j = 3:6
                if leg_list{i}(j) == 0
                    leg_list{i}(j) = -external_counter;
                    external_counter = external_counter + 1;
                end
            end
        end
    end

    map.external_legs = external_counter - 1;
    map.leg_list = leg_list;

    map.map = "true";

    if internal == 0
        if opts.boundary_matrix == 1
            map2 = PEPO.create_map(map.num_map(1:end - opts.v_cyclic, 1:end - opts.h_cyclic), opts, 1);

            boundary_map = map;
            boundary_map.boundary_map = 1;

            map = map2;
            map.boundary_map = 0;

            return;
        else
            boundary_map = map;
            return; %normal non cyclic map
        end
    else
        return;
    end
end
