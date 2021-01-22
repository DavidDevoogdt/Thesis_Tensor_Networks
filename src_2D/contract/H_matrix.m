function H_matrix = H_matrix(obj, map)

    d = obj.dim;
    H_1 = obj.H_1_tensor;
    H_2 = obj.H_2_tensor; 

    H = zeros(dimension_vector(d, 2 * map.N));

    Itensor = reshape(eye(d), [d, d, 1, 1, 1, 1]);

    %first do all single site contributions
    for l = 1:map.N

        tensor_list = cell(1, map.N);
        tensor_list(:) = {Itensor};
        tensor_list{l} = H_1;

        leg_list_copy = map.leg_list(1, :);

        for i = 1:map.N
            leg_list_copy{i} = leg_list_copy{i}(1:2);
        end

        H = H + ncon(tensor_list, leg_list_copy);
    end

    %do all horizontal H12

    for ii = 1:map.num_h_bonds

        arr = map.h_bonds{ii};
        n1 = arr(1);
        n2 = arr(2);

        if n1 ~= n2
            leg_list_copy = map.leg_list(1, :);

            index_list_n1 = map.leg_list{n1};
            index_list_n2 = map.leg_list{n2};

            leg_list_copy(max(n2, n1)) = []; %remove element

            %do ij
            new_list = [0, 0, 0, 0];

            new_list([1, 3]) = index_list_n1([1, 2]);
            new_list([2, 4]) = index_list_n2([1, 2]);

            for s = 1:map.N - 1
                leg_list_copy{s} = leg_list_copy{s}(1:2);
            end

            leg_list_copy{min(n1, n2)} = new_list;

            tensor_list = cell(1, map.N - 1);
            tensor_list(:) = {Itensor};
            tensor_list{min(n1, n2)} = H_2;

            H = H + ncon(tensor_list, leg_list_copy);
        end
    end

    %do all vertical H2
    for ii = 1:map.num_v_bonds

        arr = map.v_bonds{ii};
        n1 = arr(1);
        n2 = arr(2);

        if n1 ~= n2
            leg_list_copy = map.leg_list(1, :);

            index_list_n1 = map.leg_list{n1};
            index_list_n2 = map.leg_list{n2};

            leg_list_copy(max(n1, n2)) = []; %remove element

            for s = 1:map.N - 1
                leg_list_copy{s} = leg_list_copy{s}(1:2);
            end

            %do ij
            new_list = [0, 0, 0, 0];

            new_list([1, 3]) = index_list_n1([1, 2]);
            new_list([2, 4]) = index_list_n2([1, 2]);

            leg_list_copy{min(n1, n2)} = new_list;

            tensor_list = cell(1, map.N - 1);
            tensor_list(:) = {Itensor};
            tensor_list{min(n1, n2)} = H_2;

            H = H + ncon(tensor_list, leg_list_copy);
        end
    end

    H_matrix = reshape(H, [d^(map.N), d^(map.N)]);
end

