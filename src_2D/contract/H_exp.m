function [H, lnprefact] = H_exp(obj, map, lnprefact, improve_ln_pre_factor)

    %pos_map
    % array with ones and zeros where contraction mpo are
    %should be conected
    d = obj.dim;
    H_1 = obj.H_1_tensor;
    H_2 = obj.H_2_tensor; % +...
    %       0.5* reshape( ncon( {H_1,eye(d)}, {[-1,-3],[-2,-4]}), [d,d,d,d])+...
    %       0.5* reshape( ncon( {eye(d),H_1}, {[-1,-3],[-2,-4]}), [d,d,d,d]);

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

    %lnprefact =  abs(svds(H_matrix,1))/map.N;

    %H_expo = expm( H_matrix )/exp(lnprefactor*map.N  );

    if nargin < 3

        lnprefact = obj.nf;

    end

    if nargin < 4
        improve_ln_pre_factor = false;
    end

    %H_expo_2 = expm(H_matrix)/( exp(lnprefactor*map.N));

    d_log_nf = 1;

    while abs(d_log_nf) > 1e-10
        H_matrix_2 = H_matrix - eye(d^(map.N)) * map.N * lnprefact;

        H_expo = expm(H_matrix_2);
        d_log_nf = log(svds(H_expo, 1)) / map.N;

        lnprefact = lnprefact + d_log_nf;

        if ~improve_ln_pre_factor
            break
        end

    end

    lnprefact = lnprefact - d_log_nf;

    %Z= H_expo-H_expo_2;

    H = reshape(H_expo, dimension_vector(d, 2 * map.N));
end
