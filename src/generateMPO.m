classdef generateMPO
    %GENERATEMPO Summary of this class goes here
    %   Detailed explanation goes here

    properties
        dim
        H_1_tensor
        H_2_tensor
        type %construction type
        order
        testing
        nf
        MPO_cell
        left % left and right vector for matrix mpo, virtual level for cell mpo
        right
        MPO_type %internal representation, 'matrix' or 'cell'
        type_opts % opts specific for construction of certain type
        max_index
        current_max_index
        nf_correction
        inverse_MPO
        inverse_MPO_left
        inverse_MPO_right
    end

    properties (Access = private)
        H_exp_cell %store exp H for future uses
        H_exp_cell_cyclic

    end

    methods
        function obj = generateMPO(d, H_1_tensor, H_2_tensor, type, order, type_opts, opts)

            obj.dim = d;
            obj.H_1_tensor = H_1_tensor;
            obj.H_2_tensor = H_2_tensor;


            obj.type = type;

            if type ~= 0 %type 0 is used only for H_exp

                obj.type_opts = type_opts;

                obj.order = order;

                p = inputParser;
                addParameter(p, 'testing', 0)
                addParameter(p, 'invMPO', 1) %
                addParameter(p, 'MPO_type', 'matrix')
                parse(p, opts)


                obj.inverse_MPO = p.Results.invMPO;

                obj.MPO_type = p.Results.MPO_type;
                obj.testing = p.Results.testing;


                S = sum(abs(svds(H_exp(obj, 1, 1))));
                obj.nf = sqrt(S);


                obj = obj.makeMPO();
            end


            obj.H_exp_cell = cell(2, 1);
            obj.H_exp_cell_cyclic = cell(2, 1);
        end

        function obj = makeMPO(obj)
            switch obj.type
                case 0
                    %not really a type, just used to generate H_exp
                case 2
                    obj = obj.type_02();
                case 3
                    obj = obj.type_03();
                case 4
                    obj = obj.type_04();
                case 1
                    obj = obj.type_04();
                case 5
                    obj = obj.type_05();

                otherwise
                    error("unknown type")
            end
        end

        function obj = type_02(obj)

            d = obj.dim;

            orig_type = obj.MPO_type;

            obj.MPO_type = "cell"; %for creation,type iscasted later

            unnorm = reshape(expm(obj.H_1_tensor), [1, d, d, 1]);

            obj.MPO_cell = cell(obj.order, obj.order);

            obj.MPO_cell{0+1, 0+1} = unnorm / obj.nf;

            obj.current_max_index = 0;
            obj.max_index = obj.order;

            %all other orders
            for N = 1:obj.order
                construct_level(N)
            end


            for N = obj.order:-1:1
                err = reduce_dim(N);
                if err == 1
                    break;
                end
            end


            if orig_type == "matrix"
                mpo_cell_2_matrix();

            end

            function construct_level(N)
                %step 4 :
                % 0--|--1--|--2--|--3--|--4--|--0 = exp(H12+H23+h34)- (0--|--|--|--|--|--0)
                %N=4

                ldim = (1 + d^(2 * (N - 1))) / (1 + d^2);

                RHS_Tensor = H_exp(obj, N) / (obj.nf^(N + 1)) - obj.contract_mpo(N);
                RHS_Matrix_site = reshape(permute(RHS_Tensor, site_ordering_permute(N + 1)), ... .
                    [d^(2 * (N - 1)), d^4]);


                [part, ~] = invert_left(obj, RHS_Matrix_site, N-1);


                part = reshape(part, [d^(2 * N), d^2]);


                [U, S, V] = svd(part);

                if N ~= obj.order
                    left_i = U;
                    right_i = S * V';


                    obj = obj.register_left_inv_mpo(N-1, U');


                    obj.MPO_cell{N, N+1} = reshape(left_i, [d^(2 * (N - 1)), d, d, d^(2 * N)]);
                    obj.MPO_cell{N+1, 1} = reshape(right_i, [d^(2 * N), d, d, 1]);
                else %truncate to dim d^2
                    left_i = U * S;
                    right_i = V';

                    obj.MPO_cell{N, N+1} = reshape(left_i, [d^(2 * (N - 1)), d, d, d^2]);
                    obj.MPO_cell{N+1, 1} = reshape(right_i, [d^2, d, d, 1]);

                end

                obj.current_max_index = N;

            end

            function err = reduce_dim(N)

                O34 = obj.MPO_cell{N, N+1};
                O23 = obj.MPO_cell{N-1, N};
                O30 = obj.MPO_cell{N, 1};

                N_l = size(O23, 1);
                N_m = size(O23, 4);
                N_r = size(O34, 4);

                err = 0;

                q = d^2 * min(N_r+1, N_l);

                if N_m <= q
                    err = 1;
                    return
                end


                P1 = reshape(ncon({O23, O34}, {[-1, -2, -3, 1], [1, -4, -5, -6]}), [N_l * d^2, d^2 * N_r]);
                P2 = reshape(ncon({O23, O30}, {[-1, -2, -3, 1], [1, -4, -5, -6]}), [N_l * d^2, d^2]);


                [U, V, X, C, S] = gsvd(P1', P2');

                %undetermined dimension is lower than original
                O23_2 = reshape(X, N_l, d, d, q);
                O34_2 = reshape((U*C)', q, d, d, N_r);
                O30_2 = reshape((V*S)', q, d, d, 1);


                obj.MPO_cell{N, N+1} = O34_2;
                obj.MPO_cell{N-1, N} = O23_2;
                obj.MPO_cell{N, 1} = O30_2;

            end

            function [T, totaldimension] = mpo_cell_2_matrix()

                size_arr = zeros(obj.order+1, 1);

                for i = 1:obj.order + 1
                    size_arr(i) = size(obj.MPO_cell{i, 1}, 1);
                end

                start_index = zeros(obj.order+2, 1);
                start_index(1) = 1;
                ind = 1;
                for i = 2:obj.order + 2
                    ind = ind + size_arr(i-1);
                    start_index(i) = ind;
                end


                totaldimension = start_index(end) - 1;


                T = zeros(totaldimension, d, d, totaldimension);
                %first do tensors O_N_0
                for i = 1:obj.order + 1
                    T(start_index(i):start_index(i + 1)-1, :, :, 1) = obj.MPO_cell{i, 1};
                end

                for i = 1:obj.order
                    T(start_index(i):start_index(i + 1)-1, :, :, start_index(i + 1):start_index(i + 2)-1) = obj.MPO_cell{i, i+1};
                end

                %K=reshape(T(:,1,1,:),[45,45]);


                obj.left = zeros(1, totaldimension);
                obj.left(1) = 1;
                obj.right = zeros(totaldimension, 1);
                obj.right(1) = 1;

                obj.MPO_type = 'matrix';

                obj.MPO_cell = T;

            end

        end

        function obj = type_04(obj)
            % this type generates n--|--n blocks during the expansion
            %order n means explicit calculation op to order free bonds
            d = obj.dim;

            p = inputParser;
            addParameter(p, 'method', "diag")
            addParameter(p, 'single_threshold', 1e-12)
            addParameter(p, 'max_err', 1e0)
            addParameter(p, 'to_matrix', 0)

            parse(p, obj.type_opts)


            if mod(obj.order, 2) == 1
                obj.max_index = (obj.order + 1) / 2;
            else
                obj.max_index = obj.order / 2;
            end
            obj.current_max_index = 0;


            obj.MPO_cell = cell(obj.max_index+1, obj.max_index+1);
            obj.MPO_type = "cell";

            if obj.inverse_MPO %todo get real bounds
                obj.inverse_MPO_left = cell(obj.max_index+1, 1);
                obj.inverse_MPO_right = cell(obj.max_index+1, 1);
            end


            obj.MPO_cell{1, 1} = reshape(expm(obj.H_1_tensor)/obj.nf, [1, d, d, 1]);


            %N=number of free bonds
            for N = 1:obj.order
                if mod(N, 2) == 1
                    obj.current_max_index = (N - 1) / 2; %used to contract the tensor
                    err = double_update(N);
                    if err == 1
                        break;
                    end
                else

                    obj.current_max_index = N / 2;
                    err = single_update(N);
                    if err == 1
                        break;
                    end

                end
            end


            obj.max_index = obj.current_max_index;

            if p.Results.to_matrix == 1
                [obj.MPO_cell, total_dim] = mpo_cell_2_matrix();
                obj.MPO_type = "matrix";
            end

            function err = double_update(N)


                %step similar to
                % 0--|--1--|--2--|--1--|--0 = exp(H_12+H_23+H_34) - (0--|--|--|--|--0)


                RHS_Tensor = H_exp(obj, N) / (obj.nf^(N + 1)) - obj.contract_mpo(N);


                RHS_Matrix = reshape(permute(RHS_Tensor, site_ordering_permute(N + 1)), ...
                    dimension_vector(d^2, N + 1)); %group per ij index

                M = obj.current_max_index;


                res = reshape(RHS_Matrix, [d^(2 * M), d^4 * d^(2 * M)]);

                [x, left_dim] = invert_left(obj, res, M);
                x = reshape(x, [left_dim * d^4, d^(2 * M)]);
                [y, right_dim] = invert_right(obj, x, M);

                new_parts = reshape(y, [left_dim, d^2, d^2, right_dim]);

                new_parts_sym = reshape(permute(new_parts, [1, 2, 4, 3]), [left_dim * d^2, d^2 * right_dim]);


                %%%%%%%%%%%%
                switch p.Results.method 
                    case "svd"
              
                        [U, S, V] = svd(new_parts_sym);
                        
                                        
                        index_arr = diag(S) > p.Results.single_threshold;
                        middle_dim = sum(index_arr);

                        U2 = U(:, 1:middle_dim);
                        S2 = S(1:middle_dim, 1:middle_dim);
                        V2 = V(:, 1:middle_dim);
                        

                        V2_dag_perm = reshape(permute(reshape(V2', [middle_dim, right_dim, d, d]), [1, 3, 4, 2]), [ middle_dim,right_dim*d^2]);
                        V2_perm = reshape(permute(reshape(V2, [right_dim, d, d, middle_dim]), [2, 3, 1, 4]), [right_dim*d^2, middle_dim]);


                        sqrt_S = diag(diag(S2).^0.5);
                        sqrt_S_inv = diag(diag(S2).^(-0.5));
                        
                        U_inv = sqrt_S_inv*U2';
                        V_inv = V2_perm*sqrt_S_inv;
                            
                        
                        U = U2 * sqrt_S;
                        V = sqrt_S * V2_dag_perm;
                        
                    case "diag"
                        [U, S] = eig(new_parts_sym);
                        middle_dim = size(S,1);
                        U_inv= U^-1;  %expensive
                        
                        V2_dag_perm = reshape(permute(reshape(U_inv, [middle_dim, right_dim, d, d]), [1, 3, 4, 2]), [middle_dim, middle_dim]);
                        V2_perm = reshape(permute(reshape(U, [right_dim, d, d, middle_dim]), [2, 3, 1, 4]), [middle_dim, middle_dim]);

                        sqrt_S_vect = diag(S).^0.5;
                        sqrt_S = diag(sqrt_S_vect);
                        sqrt_S_inv = diag(sqrt_S_vect.^(-1));
                        
                        U_inv =  sqrt_S_inv*U_inv;
                        V_inv = V2_perm*sqrt_S_inv;
                        
                        U = U*sqrt_S;
                        V = sqrt_S*V2_dag_perm;
                    otherwise
                        error("unknown method")
                        
                end


                %%%%%%%%%%%%%%%%%%%


                O_l = reshape(U, [left_dim, d, d, middle_dim]);
                %O_r = reshape(right_j, [middle_dim, d, d,right_dim]);
                O_r = reshape(V, [middle_dim, d, d, right_dim]);

                obj = obj.register_left_inv_mpo(M, U_inv);
                obj = obj.register_right_inv_mpo(M, V_inv);


                hexp2 = H_exp(obj, N, 0, 1); %cyclic error
                hexp = hexp2 / (obj.nf^(N + 1));

                %nf_corr = sum(  abs(svds(reshape(hexp,[d^(2*M+2),d^(2*M+2)]) ))).^(1/(N+2));


                RHS_Tensor_2 = hexp - obj.contract_mpo(N, 0, 1);


                RHS_Matrix_2 = reshape(permute(RHS_Tensor_2, site_ordering_permute(N + 1)), ...
                    dimension_vector(d^2, N + 1));


                obj.MPO_cell{M+1, M+2} = O_l;
                obj.MPO_cell{M+2, M+1} = O_r;
                obj.current_max_index = M + 1;

                RHS_Tensor_3 = hexp - obj.contract_mpo(N, 0, 1);
                RHS_Matrix_3 = reshape(permute(RHS_Tensor_3, site_ordering_permute(N + 1)), ...
                    dimension_vector(d^2, N + 1));


                RHS_2_reshaped = reshape(RHS_Matrix_2, [d^(2 * M + 2), d^(2 * M + 2)]);
                RHS_3_reshaped = reshape(RHS_Matrix_3, [d^(2 * M + 2), d^(2 * M + 2)]);

                cond1 = svds(RHS_2_reshaped, 6);
                cond2 = svds(RHS_3_reshaped, 6);

                p1 = 2;
                sum1 = (sum(cond1.^p1))^(1 / p1);
                sum2 = (sum(cond2.^p1))^(1 / p1);

                condres = sum2 / sum1;

                err = 0;

                %fprintf("-%.4e %.4e-",sumcond,maxcond)

                if condres > p.Results.max_err
                    fprintf("-%.4e -", condres)
                    obj.MPO_cell{M+1, M+2} = {};
                    obj.MPO_cell{M+2, M+1} = {};
                    obj.current_max_index = M;
                    fprintf("stop %d", N);
                    err = 1;
                end

            end

            function err = single_update(N)


                RHS_Tensor = H_exp(obj, N) / (obj.nf^(N + 1)) ...
                    -obj.contract_mpo(N);

                RHS_Matrix = reshape(permute(RHS_Tensor, site_ordering_permute(N + 1)), ...
                    dimension_vector(d^2, N + 1)); %group per ij index


                %search x st
                % left*x = res
                % y*right = x
                M = obj.current_max_index;


                res = reshape(RHS_Matrix, [d^(2 * M), d^(2 * M + 2)]);

                [x, left_dim] = invert_left(obj, res, M);
                x = reshape(x, [left_dim * d^2, d^(2 * M)]);
                [y, right_dim] = invert_right(obj, x, M);


                err = 0;

                %cyclic error reduction

                hexp2 = H_exp(obj, N+1, 0, 1); %cyclic error
                hexp = hexp2 / (obj.nf^(N + 2));

                %nf_corr = sum(  abs(svds(reshape(hexp,[d^(2*M+2),d^(2*M+2)]) ))).^(1/(N+2));


                RHS_Tensor_2 = hexp - obj.contract_mpo(N+1, 0, 1);


                RHS_Matrix_2 = reshape(permute(RHS_Tensor_2, site_ordering_permute(N + 2)), ...
                    dimension_vector(d^2, N + 2));


                obj.MPO_cell{M+1, M+1} = reshape(y, [left_dim, d, d, right_dim]);

                RHS_Tensor_3 = hexp - obj.contract_mpo(N+1, 0, 1);
                RHS_Matrix_3 = reshape(permute(RHS_Tensor_3, site_ordering_permute(N + 2)), ...
                    dimension_vector(d^2, N + 2));


                RHS_2_reshaped = reshape(RHS_Matrix_2, [d^(2 * M + 2), d^(2 * M + 2)]);
                RHS_3_reshaped = reshape(RHS_Matrix_3, [d^(2 * M + 2), d^(2 * M + 2)]);


                cond1 = svds(RHS_2_reshaped, 10);
                cond2 = svds(RHS_3_reshaped, 10);

                p1 = 2;
                sum1 = (sum(cond1.^p1))^(1 / p1);
                sum2 = (sum(cond2.^p1))^(1 / p1);
                condres = sum2 / sum1;


                %fprintf("-%.4e %.4e-",sumcond,maxcond)

                if condres > p.Results.max_err

                    fprintf("-%.4e -", condres)
                    obj.MPO_cell{M+1, M+1} = {};
                    fprintf("stop %d", N);
                    err = 1;
                end

            end

            function [T, totaldimension] = mpo_cell_2_matrix()

                size_arr = zeros(obj.max_index+1, 1);

                size_arr(1) = 1;

                for i = 1:obj.max_index
                    size_arr(i+1) = size(obj.MPO_cell{i, i + 1}, 4);
                end

                start_index = zeros(obj.max_index+2, 1);
                start_index(1) = 1;
                ind = 1;
                for i = 2:obj.max_index + 2
                    ind = ind + size_arr(i-1);
                    start_index(i) = ind;
                end


                totaldimension = start_index(end) - 1;


                T = zeros(totaldimension, d, d, totaldimension);
                %first do tensors O_N_0
                for i = 1:obj.max_index + 1
                    if size(obj.MPO_cell{i, i}, 1) ~= 0
                        T(start_index(i):start_index(i + 1)-1, :, :, start_index(i):start_index(i + 1)-1) = obj.MPO_cell{i, i};
                    end


                    if i > 1
                        T(start_index(i - 1):start_index(i)-1, :, :, start_index(i):start_index(i + 1)-1) = obj.MPO_cell{i-1, i};
                        T(start_index(i):start_index(i + 1)-1, :, :, start_index(i - 1):start_index(i)-1) = obj.MPO_cell{i, i-1};
                    end
                end


                %K=reshape(T(:,1,1,:),[45,45]);


                obj.left = zeros(1, totaldimension);
                obj.left(1) = 1;
                obj.right = zeros(totaldimension, 1);
                obj.right(1) = 1;

                obj.MPO_type = 'matrix';

                obj.MPO_cell = T;

            end


        end

        function obj = type_03(obj)
            % Make a cell from O. It holds the tensor elements
            % every entry is 4d nxdxdxm with n and m the bond dimension for the
            % corresponding bond. dimension x1 at the end not shown by matlab
            %this type generates no 1--|--1 and 2--|--2 blocks
            %uses virtual levels:
            %0--|--1--|--0
            %0--|--1'--|--2'--|--0
            %0--|--1''--|--2''--|--3''--|--0


            %representation blocks in block matrix with dims

            %
            %       1         d^2     d^2    d^2      d^2     d^4         d^2
            %
            %  1     00        01'    |01''   0       |01'''   0           0       |
            %  d^2   1'0       0      |0      0       |0       0           0       |
            %       __________________|               |                            |
            %  d^2   0         0       0       1''2'' |0       0           0       |
            %  d^2   2''0      0       0       0      |0       0           0       |
            %       __________________________________|0       0           0       |
            %  d^2   0         0       0       0       0       1'''2'''    0       |
            %  d^4   0         0       0       0       0       0           2'''3'''|
            %  d^2   3'''0     0       0       0       0       0           0       |
            %       _______________________________________________________________|


            d = obj.dim;

            total_dim = 1;

            for k = 1:obj.order
                for i = 1:k
                    total_dim = total_dim + internal_dim(i, k, d);
                end
            end


            %total_dim = 1+ d^2 + (d^2+d^2)+ (d^2+d^4+d^2)+(d^2+d^4+d^4+d^2);


            p = inputParser;
            addParameter(p, 'SparseArr', 1)

            parse(p, obj.type_opts)

            if p.Results.SparseArr == 1
                sparsem = ndSparse(sparse(total_dim, d^2 * total_dim));
                obj.MPO_cell = reshape(sparsem, [total_dim, d, d, total_dim]);
            else
                obj.MPO_cell = zeros(total_dim, d, d, total_dim);
            end


            left_vect = zeros(1, total_dim);
            left_vect(1) = 1;
            right_vect = zeros(total_dim, 1);
            right_vect(1) = 1;


            %obj.MPO_cell = zeros(total_dim,d,d,total_dim);
            obj.MPO_type = "matrix";
            obj.left = left_vect;
            obj.right = right_vect;


            % 0
            unnorm = expm(obj.H_1_tensor);
            %obj.nf = trace(unnorm);

            O_00_0 = reshape(unnorm/obj.nf, [1, d, d, 1]);

            obj.MPO_cell = add_block_to_tensor(obj.MPO_cell, 0, 0, 0, O_00_0, d);

            %obj.MPO_cell{0,0} = O_00_0;

            %step 1:
            % 0--|--1'--|--0 = exp(H_12) - (0--|--|--0 )
            N = 1; %number accents = number of free bonds


            RHS_Tensor = H_exp(obj, N) / (obj.nf^(N + 1)) - obj.contract_mpo(N);
            RHS_Matrix_site = reshape(permute(RHS_Tensor, site_ordering_permute(N + 1)), ...
                [d^2, d^2]); %ready to svd

            O_01_1 = reshape(eye(d^2), [1, d, d, d^2]);
            O_10_1 = reshape(RHS_Matrix_site, [d^2, d, d, 1]);

            obj.MPO_cell = add_block_to_tensor(obj.MPO_cell, 0, 1, N, O_01_1, d);
            obj.MPO_cell = add_block_to_tensor(obj.MPO_cell, 1, 0, N, O_10_1, d);


            if obj.testing == 1
                err = tensor_norm(ncon({O_01_1, O_10_1}, {[-1, -2, -4, 1], [1, -3, -5, -6]})-RHS_Tensor);
                fprintf("err 01 = %d\n", err);
            end

            %all other orders
            for N = 2:obj.order
                construct_level(N)
            end


            function construct_level(N)
                %step 4 :
                % 0--|--1''''--|--2''''--|--3''''--|--4''''--|--0 = exp(H12+H23+h34)- (0--|--|--|--|--|--0)
                %N=4

                %N=3
                %0--|--1'''--|--2'''--|--3'''---|---0

                if mod(N, 2) == 0
                    M_l = N / 2;
                    M_r = N / 2;
                else
                    M_l = (N - 1) / 2;
                    M_r = (N - 1) / 2 + 1;
                end

                RHS_Tensor = H_exp(obj, N) / (obj.nf^(N + 1)) - obj.contract_mpo(N);
                RHS_Matrix_site = reshape(permute(RHS_Tensor, site_ordering_permute(N + 1)), ...
                    [d^(2 * M_l), d, d, d^(2 * M_r)]);

                for k = 0:M_l - 1
                    obj.MPO_cell = add_block_to_tensor(obj.MPO_cell, k, k+1, N, reshape(eye(d^(2 * (k + 1))), [d^(2 * k), d, d, d^(2 * (k + 1))]), d);
                end

                obj.MPO_cell = add_block_to_tensor(obj.MPO_cell, M_l, M_l+1, N, reshape(RHS_Matrix_site, [d^(2 * M_l), d, d, d^(2 * M_r)]), d);

                for k = M_l + 1:M_l + M_r - 1
                    obj.MPO_cell = add_block_to_tensor(obj.MPO_cell, k, k+1, N, reshape(eye(d^(2 * (N - k + 1))), [d^(2 * (N - k + 1)), d, d, d^(2 * (N - k))]), d);
                end

                obj.MPO_cell = add_block_to_tensor(obj.MPO_cell, N, 0, N, reshape(eye(d^2), [d^2, d, d, 1]), d);
            end

            function y = internal_dim(i, k, d)
                y = d^(2 * min(i, k - i + 1));
            end

            %test function for type_03
            %             function test_add_block_to_tensor
            %                 d=2;
            %                 dim = 1+d^2+(d^2+d^2) + (d^2+d^4+d^2);
            %                 O = zeros(dim,d,d,dim);
            %
            %                 b_00_1 = ones(1,d,d,1)*0.5;
            %                 O = add_block_to_tensor(O,0,0,0,b_00_1,d);
            %
            %                 b_01_1 = ones(1,d,d,d^2);
            %                 O = add_block_to_tensor(O,0,1,1,b_01_1,d);
            %
            %                 b_10_1 = ones(d^2,d,d,1)*2;
            %                 O = add_block_to_tensor(O,1,0,1,b_10_1,d);
            %
            %                 b_01_2 = ones(1,d,d,d^2)*3;
            %                 O = add_block_to_tensor(O,0,1,2,b_01_2,d);
            %
            %                 b_12_2 = ones(d^2,d,d,d^2)*5;
            %                 O = add_block_to_tensor(O,1,2,2,b_12_2,d);
            %
            %                 b_20_2 = ones(d^2,d,d,1)*4;
            %                 O = add_block_to_tensor(O,2,0,2,b_20_2,d);
            %
            %
            %                 b_01_3 = ones(1,d,d,d^2)*6;
            %                 O = add_block_to_tensor(O,0,1,3,b_01_3,d);
            %
            %                 b_12_3 = ones(d^2,d,d,d^4)*7;
            %                 O = add_block_to_tensor(O,1,2,3,b_12_3,d);
            %
            %                 b_23_3 = ones(d^4,d,d,d^2)*8;
            %                 O = add_block_to_tensor(O,2,3,3,b_23_3,d);
            %
            %                 b_30_3 = ones(d^2,d,d,1)*9;
            %                 O = add_block_to_tensor(O,3,0,3,b_30_3,d);
            %
            %                 Z= reshape(O(:,1,1,:),[dim,dim]);
            %
            %             end
            %

            %representation blocks in block matrix.

            %
            %        1         d^2     d^2    d^2      d^2     d^4         d^2
            %
            %  1     00        01'    |01''   0       |01'''   0           0       |
            %  d^2   1'0       0      |0      0       |0       0           0       |
            %       __________________|               |                            |
            %  d^2   0         0       0       1''2'' |0       0           0       |
            %  d^2   2''0      0       0       0      |0       0           0       |
            %       __________________________________|0       0           0       |
            %  d^2   0         0       0       0       0       1'''2'''    0       |
            %  d^4   0         0       0       0       0       0           2'''3'''|
            %  d^2   3'''0     0       0       0       0       0           0       |
            %       _______________________________________________________________|


            %add block  i(k)--|--j(k) to the matrix
            %O (dim,d,d,dim)
            function O = add_block_to_tensor(O, i, j, k, block, d)

                block_start = 1; %00 block

                if i == 0 && j == 0
                    O(1, :, :, 1) = block(:, :, :, :);
                    return;
                end

                %todo make explicit formula
                for s = 1:k - 1
                    for t = 1:s
                        block_start = block_start + internal_dim(t, s, d);
                    end
                end

                if i == 0
                    y = block_start + 1;
                    x = 1;

                else
                    if j == 0
                        x = block_start + 1;
                        for s = 1:i - 1
                            x = x + internal_dim(s, k, d);
                        end
                        y = 1;
                    else
                        x = block_start + 1;
                        for s = 1:i - 1
                            x = x + internal_dim(s, k, d);
                        end

                        y = block_start + 1;
                        for t = 1:j - 1
                            y = y + internal_dim(t, k, d);
                        end
                    end
                end

                xdim = size(block, 1);
                ydim = size(block, 4);
                O(x:x+xdim-1, :, :, y:y+ydim-1) = O(x:x+xdim-1, :, :, y:y+ydim-1) + block(:, :, :, :);

                % k is subspace dim
                function y = internal_dim(i, k, d)
                    y = d^(2 * min(i, k - i + 1));
                end
            end
        end

        function obj = type_05(obj)
            %todo matrix not reallly necessary and large dimension, make
            %normal cell mpo

            %this mpo limits the possibilities and works with unitary transformations
            %
            %  unitary: (S_(n n+1)^(i j alpha j)_beta) ' = T_(n+1 n)^ alpha _ (i j beta )
            % matrix is constructed st only these combinations can happen
            %  S01--S12--...--S(n-1 n)--D_nn--n--S(n-1 n)'--...--S21'--S10'    with D a real diagonal matrix
            %  S01--S12--...--S(n-1 n)--(T_nn--)^m--S(n-1 n)'--...--S21'--S10' with
            %  T_nn a tesnsor

            % cosntruction matrix,  ' means hermitian conjugate. Example for blocks uo
            % until 2, can be generalised for any n
            %      B 0 |   B1     |     B2     |     B3    |    B4

            % dims
            %        1   d^2  d^4   d^2   d^4    d^2  d^4    d^2      d^4
            %
            % 1    T_00| S01  0   |-2*S01  0   | S01  0    | S01*D11   0
            %      ____|          |            | 0    S12  | 0         0
            % d^2  S01'  0    S12 | 0 -2*S12   | 0    0    | 0         0
            % d^4  0     S12' 0   | 0      0   | 0    0    | 0         0
            %      _______________| 0      0   | 0    0    | 0         0
            % d^2  S01'  0    0     0      0   | 0    0    | 0         0
            % d^4  0     S12' 0     0      0   | 0    0    | 0         0
            %      ____________________________| 0    0    | 0         0
            % d^2  S01'  0    0     0      0     T11  0    | 0         0
            % d^4  0     S12' 0     0      0     0    T22  | 0         0
            %      ________________________________________| 0         0
            % d^2  S01'  0    0     0      0     0    0      0         (1/D11)*S12*D22
            % d^4  0     0    0     0      0     0    0      S12'      0

            %intuition matrix: B4: create blocks with D,  B3 with B1: create T_nn
            %sequences. B1 is is used to create S(n n+1) sequence beforen T_nn^m
            % block B2 corrects for spurious multiplications in B1 and B3 that do not
            % involve any T_nn


            obj.nf = 3 * obj.nf;

            p = inputParser;
            addParameter(p, 'method', "diag")
            addParameter(p, 'SparseArr', 1)

            parse(p, obj.type_opts)


            %setup
            d = obj.dim;

            if mod(obj.order, 2) == 1
                obj.max_index = (obj.order + 1) / 2;
            else
                obj.max_index = obj.order / 2;
            end

            obj.current_max_index = 0;

            total_dim = 1 + 4 * (((d^2)^(obj.max_index + 1) - 1) / (d^2 - 1) - 1);

            obj.left = zeros(1, total_dim);
            obj.left(1) = 1;
            obj.right = zeros(total_dim, 1);
            obj.right(1) = 1;

            if p.Results.SparseArr == 1
                sparsem = ndSparse(sparse(total_dim, d^2 * total_dim));
                obj.MPO_cell = reshape(sparsem, [total_dim, d, d, total_dim]);
            else
                obj.MPO_cell = zeros(total_dim, d, d, total_dim);
            end


            %

            %U_cell = cell(1, obj.max_index);
            %V_cell = cell(1, obj.max_index);

            sqrt_Dn_l = [1];
            sqrt_Dn_r = [1];

            % 00 block
            N = 0;


            obj.MPO_type = "matrix";

            T00 = reshape(expm(obj.H_1_tensor)/obj.nf, [1, d, d, 1]);
            add_block_Tn(N, T00, obj.order, d);

            %other blocks
            %N=number of free bonds
            for N = 1:obj.order
                if mod(N, 2) == 1
                    %obj.current_max_index = (N-1)/2;  %used to contract the tensor
                    [sqrt_Dn_l, sqrt_Dn_r] = double_update(N, sqrt_Dn_l, sqrt_Dn_r);
                else
                    %obj.current_max_index = N/2;
                    single_update(N)
                end

            end


            %n current num of bonds
            %N max free bonds (order)
            function [sqrt_Dn_l_inv, sqrt_Dn_r_inv] = double_update(N, sqrt_Dnm_l_inv, sqrt_Dnm_r_inv)
                %  S01--S12--...--S(n-1 n)--D_nn--n--S(n-1 n)'--...--S21'--S10' with D a real diagonal matrix


                RHS_Tensor = H_exp(obj, N) / (obj.nf^(N + 1)) - obj.contract_mpo(N);


                RHS_Matrix = reshape(permute(RHS_Tensor, site_ordering_permute(N + 1)), ...
                    dimension_vector(d^2, N + 1)); %group per ij index

                M = obj.current_max_index;


                res = reshape(RHS_Matrix, [d^(2 * M), d^4 * d^(2 * M)]);

                %[x, left_dim] = invert_left(obj, res, M, U_cell);
                [x, left_dim] = invert_left(obj, res, M);
                x = reshape(x, [left_dim * d^4, d^(2 * M)]);
                %[y, right_dim] = invert_right(obj, x, M, V_cell);
                [y, right_dim] = invert_right(obj, x, M);


                new_parts = reshape(y, [left_dim, d^2, d^2, right_dim]);


                new_parts_sym = reshape(permute(new_parts, [1, 2, 4, 3]), [left_dim * d^2, d^2 * right_dim]);


                %new_parts= reshape(  new_parts , [left_dim*d^2,d^2*right_dim]);

                switch p.Results.method 
                    case "svd"
              
                        [U2, S, V2] = svd(new_parts_sym);

                        V2_dag = V2';
                        V2_inv = V2;
                        U2_inv = U2';
                        
                       
                        
                    case "diag"
                        [U2, S] = eig(new_parts_sym);
                        
                        V2_dag = U2^-1;
                        V2_inv = U2;
                        U2_inv = V2_dag;
                    otherwise
                        error("unknown method")
                        
                end
                
             
                 
                        middle_dim = size(S,1);
                        

                        V2_dag_perm = reshape(permute(reshape(V2_dag, [middle_dim, right_dim, d, d]), [1, 3, 4, 2]), [middle_dim, middle_dim]);
                        V2_perm = reshape(permute(reshape(V2_inv, [right_dim, d, d, middle_dim]), [2, 3, 1, 4]), [middle_dim, middle_dim]);

                        U_inv = U2_inv;
                        V_inv = V2_perm;
                            
                        U = U2 ;
                        V_dag = V2_dag_perm;
                
              
                        %[U, Dn, V] = svd(new_parts);

%                         middle_dim = size(Dn, 1);
% 
%                         eps = 1e-15;
% 
                         maxval = sum(abs(diag(S)));
                         sqrtmaxval = maxval^0.5;
% 
                         diag_Dn = diag(S/maxval);

                         diag_Dn_inv = diag_Dn.^-1;

%                         diag_Dn_inv = diag_Dn;
%                         diag_Dn_inv(~mask) = diag_Dn_inv(~mask).^(-1);
% 
                        diag_vect = diag_Dn.^(1 / 2);

                         sqrt_Dn_l = diag(diag_vect);
                         sqrt_Dn_r = sqrt_Dn_l;
                         sqrt_Dn_l_inv = diag(diag_vect.^(-1));
                         sqrt_Dn_r_inv = sqrt_Dn_l_inv;


                        obj = register_left_inv_mpo(obj, M, U_inv / sqrtmaxval);
                        obj = register_right_inv_mpo(obj, M, V_inv / sqrtmaxval);

                        Sn = reshape(U*sqrtmaxval, [left_dim, d, d, middle_dim]);
                        SnD = reshape(V_dag*sqrtmaxval, [middle_dim, d, d, right_dim]);
                        %SnD =  reshape(V'*sqrtmaxval, [d^(2*M+2),d,d,d^(2*M)] ) ;
                        %U_cell{1+obj.current_max_index} = Sn;
                        %V_cell{1+obj.current_max_index} = SnD;



                add_block_05(obj.current_max_index+1, Sn, SnD, sqrt_Dn_l, sqrt_Dn_r, sqrt_Dnm_l_inv, sqrt_Dnm_r_inv, d);

                obj.current_max_index = obj.current_max_index + 1;

            end


            function single_update(N)
                %  S01--S12--...--S(n-1 n)--(T_nn--)^m--S(n n-1)'--...--S21'--S10' with


                RHS_Tensor = H_exp(obj, N) / (obj.nf^(N + 1)) ...
                    -obj.contract_mpo(N);

                RHS_Matrix = reshape(permute(RHS_Tensor, site_ordering_permute(N + 1)), ...
                    dimension_vector(d^2, N + 1)); %group per ij index

                M = obj.current_max_index;

                res = reshape(RHS_Matrix, [d^(2 * M), d^2 * d^(2 * M)]);

                
                %[x, left_dim] = invert_left(obj, res, M, U_cell);
                [x, left_dim] = invert_left(obj, res, M);
                x = reshape(x, [left_dim * d^2, d^(2 * M)]);
                %[y, right_dim] = invert_right(obj, x, M, V_cell);
                [y, right_dim] = invert_right(obj, x, M);

                new_parts = reshape(y, [left_dim, d, d, right_dim]);

                %add_block_Tn(MPO,n,Tn,max_index,d )
                add_block_Tn(obj.current_max_index, new_parts, obj.max_index, d);

            end

            %add blocks to mpo like in discription
            function add_block_05(n, Sn, SnD, sqrt_Dn_l, sqrt_Dn_r, sqrt_Dnm_l_inv, sqrt_Dnm_r_inv, d)

                %horizontal
                block_start_x = 1;

                %B=1 case
                B = 1;
                block_start_y = get_B_start(B, d, obj.max_index);
                internal_diag = geom_sum(n-1, d);
                internal_diag_m = geom_sum(n-2, d) + 1;
                obj.MPO_cell(1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2 * n) - 1), :, :, block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2 * n - 2) - 1)) = SnD;
                %mirror case for
                obj.MPO_cell(block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2 * n - 2) - 1), :, :, 1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2 * n) - 1)) = Sn;

                %                 if n==1
                %B=2 case
                B = 2;
                block_start_y = get_B_start(B, d, obj.max_index);
                internal_diag = geom_sum(n-1, d);
                internal_diag_m = geom_sum(n-2, d) + 1;
                obj.MPO_cell(1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2 * n) - 1), :, :, block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2 * n - 2) - 1)) = SnD;
                %mirror case for
                obj.MPO_cell(block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2 * n - 2) - 1), :, :, 1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2 * n) - 1)) = -2 * Sn;

                %B=3 case
                B = 3;
                block_start_y = get_B_start(B, d, obj.max_index);
                internal_diag = geom_sum(n-1, d);
                internal_diag_m = geom_sum(n-2, d) + 1;
                obj.MPO_cell(1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2 * n) - 1), :, :, block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2 * n - 2) - 1)) = SnD;
                %mirror case for
                obj.MPO_cell(block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2 * n - 2) - 1), :, :, 1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2 * n) - 1)) = Sn;


                if n == 1

                    SD = ncon({Sn, sqrt_Dn_l}, {[-1, -2, -3, 1], [1, -4]});

                    SDn = ncon({sqrt_Dn_r, SnD}, {[-1, 1], [1, -2, -3, -4]});

                    block_start_y = get_B_start(4, d, obj.max_index);
                    %MPO
                    obj.MPO_cell(1, :, :, 1+block_start_y:1+block_start_y+(d^2 - 1)) = SD;
                    obj.MPO_cell(1+block_start_y:1+block_start_y+(d^2 - 1), :, :, 1) = SDn;

                else %todo
                    block_start_y = get_B_start(4, d, obj.max_index);
                    internal_diag = geom_sum(n-1, d);
                    internal_diag_m = geom_sum(n-2, d) + 1;


                    %                     sqrt_Dnm_l_inv = sqrt_Dnm_l_inv^-1; %is diagonal
                    %                     sqrt_Dnm_r_inv = sqrt_Dnm_r_inv^-1; %is diagonal
                    %
                    SD = ncon({sqrt_Dnm_l_inv, Sn, sqrt_Dn_l}, {[-1, 1], [1, -2, -3, 2], [2, -4]});

                    SDn = ncon({sqrt_Dn_r, SnD, sqrt_Dnm_r_inv}, {[-1, 1], [1, -2, -3, 2], [2, -4]});

                    obj.MPO_cell(block_start_y+internal_diag_m:block_start_y+internal_diag_m+(d^(2 * n - 2) - 1), :, :, 1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2 * n) - 1)) = SD;
                    obj.MPO_cell(1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2 * n) - 1), :, :, block_start_y+internal_diag_m:block_start_y+internal_diag_m+(d^(2 * n - 2) - 1)) = SDn;
                end

                function y = get_B_start(B, d, N)
                    block_dim = geom_sum(N, d);

                    y = 1 + (B - 1) * block_dim;
                end

                function y = geom_sum(N, d)
                    y = ((d^2)^(N + 1) - 1) / (d^2 - 1) - 1;
                end
            end


            function add_block_Tn(n, Tn, max_index, d)

                if n == 0
                    obj.MPO_cell(1, :, :, 1) = Tn;
                    return
                end

                B = 3;
                block_start_y = get_B_start(B, d, max_index);
                internal_diag = geom_sum(n-1, d);

                z = block_start_y + internal_diag;
                obj.MPO_cell(1+z:1+z+(d^(2 * n) - 1), :, :, 1+z:1+z+(d^(2 * n) - 1)) = Tn;

                function y = get_B_start(B, d, N)
                    block_dim = geom_sum(N, d);

                    y = 1 + (B - 1) * block_dim;
                end

                function y = geom_sum(N, d)
                    y = ((d^2)^(N + 1) - 1) / (d^2 - 1) - 1;
                end
            end

            %test code for previous function, everything works fine
            %old version where t also was added in same fn
            function test_add_block
                d = 2;
                N = 3;
                total_dim = 1 + 4 * (((d^2)^(N + 1) - 1) / (d^2 - 1) - 1);
                MPO = zeros(total_dim, d, d, total_dim);

                %     add_block_05(MPO,n,Sn,SnD,Tn,            Dn,Dnm,N,d )
                MPO = add_block_05(MPO, 0, 0, 0, 0.3*ones(1, d, d, 1), 0, 0, N, 2);

                %S1
                N = 1;
                S01 = ones(1, d, d, d^2) * 0.5454;
                S10 = ones(d^2, d, d, 1) * 0.5454;
                T11 = ones(d^2, d, d, d^2) * 3.3333;
                D11 = eye(d^2);
                %     add_block_05(MPO,n,Sn,SnD,Tn,Dn,Dnm,N,d )
                MPO = add_block_05(MPO, N, S01, S10, T11, D11, 1, N, 2);

                %S2
                N = 2;
                S12 = ones(d^2, d, d, d^4) * 0.156165;
                S21 = ones(d^4, d, d, d^2) * 0.156165;
                T22 = ones(d^4, d, d, d^4) * 0.1686412;
                D22 = eye(d^4);
                %     add_block_05(MPO,n,Sn,SnD,Tn,Dn,Dnm,N,d )
                MPO = add_block_05(MPO, N, S12, S21, T22, D22, D11, N, 2);


                %S3
                N = 3;
                S23 = ones(d^4, d, d, d^6) * 0.7777;
                S32 = ones(d^6, d, d, d^4) * 0.7777;
                T33 = ones(d^6, d, d, d^6) * 0.1686412;
                D33 = eye(d^6);
                %     add_block_05(MPO,n,Sn,SnD,Tn,Dn,Dnm,N,d )
                MPO = add_block_05(MPO, N, S23, S32, T33, D33, D22, N, 2);

                %Z=reshape( MPO(:,1,1,:), [total_dim,total_dim]) ;
            end

        end

        function H_exp = H_exp(obj, N, matrix, cyclic)
            if nargin < 4
                cyclic = 0;
            end
            if nargin < 3
                matrix = 0;
            end

            in_memory = 0;

            if cyclic == 1
                if size(obj.H_exp_cell_cyclic, 1) >= N
                    if obj.H_exp_cell_cyclic{N} ~= []
                        in_memory = 1;
                    end
                end
            else
                if size(obj.H_exp_cell, 1) >= N
                    if obj.H_exp_cell{N} ~= []
                        in_memory = 1;
                    end
                end
            end

            if in_memory == 0

                % return E(H_1_2+..+H_N-1_N) in normal ordering (dimension d^N+1, basis first
                % upper legs, then lower legs
                % so this makes first the tensor T = H     x I x I
                %                                  + I x H     x I ...
                %
                %                                  + S x I x I x I
                %                                  + I x S x I x I ...
                % with S the single site operator and H the 2 site one
                % then reorders, and exponentiates
                d = obj.dim;

                H = zeros(dimension_vector(d, 2 * (N + 1)));

                %1 site operator
                for i = 1:N + 1
                    tensor_list = cell(1, N+1);
                    leg_list = cell(1, N+1);

                    for j = 1:i - 1
                        tensor_list{j} = eye(d);
                        leg_list{j} = [-j, -(N + 1 + j)];
                    end

                    tensor_list{i} = obj.H_1_tensor;
                    leg_list{i} = [-i, -(N + 1 + i)];

                    for j = i + 1:N + 1
                        tensor_list{j} = eye(d);
                        leg_list{j} = [-(j), -(N + 1 + j)];
                    end


                    H_i = ncon(tensor_list, leg_list);
                    H = H + H_i;
                end


                %2 site operator
                for i = 1:N
                    tensor_list = cell(1, N);
                    leg_list = cell(1, N);

                    for j = 1:i - 1
                        tensor_list{j} = eye(d);
                        leg_list{j} = [-j, -(N + 1 + j)];
                    end

                    tensor_list{i} = obj.H_2_tensor;
                    leg_list{i} = [-i, -(i + 1), -(N + 1 + i), -(N + i + 2)];

                    for j = i + 1:N
                        tensor_list{j} = eye(d);
                        leg_list{j} = [-(j + 1), -(N + j + 2)];
                    end

                    H_i = ncon(tensor_list, leg_list, []);
                    H = H + H_i;
                end

                %cyclic case

                if cyclic == 1
                    tensor_list = cell(1, N);
                    leg_list = cell(1, N);

                    for j = 2:N
                        tensor_list{j} = eye(d);
                        leg_list{j} = [-j, -(N + 1 + j)];
                    end


                    tensor_list{1} = obj.H_2_tensor;
                    leg_list{1} = [-(N + 1), -1, -(2 * N + 2), -(N + 2)];


                    H_i = ncon(tensor_list, leg_list, []);
                    H = H + H_i;
                end
            else
                if cyclic == 1
                    H = obj.H_exp_cell_cyclic{N};
                else
                    H = obj.H_exp_cell{N};
                end
            end


            H_exp_matrix = expm(reshape(H, [d^(N + 1), d^(N + 1)]));

            if matrix == 1
                H_exp = H_exp_matrix;
                return
            end

            if cyclic == 1
                H_exp = reshape(H_exp_matrix, dimension_vector(d, 2 * (N + 1)));
            else
                H_exp = reshape(H_exp_matrix, dimension_vector(d, 2 * (N + 1), [1, 1]));
            end

        end

        function T = contract_mpo(obj, N, matrix, cyclic, virtual_level)
            if nargin < 5
                virtual_level = 0;
            end

            if nargin < 4
                cyclic = 0;
            end
            if nargin < 3
                matrix = 0;
            end

            switch obj.MPO_type
                case "cell"

                    O = obj.MPO_cell;
                    d = obj.dim;

                    if cyclic == 1
                        total_pos = (obj.current_max_index + 1)^(N + 1) - 1;
                        T = zeros(dimension_vector(d, 2 * (N + 1)));
                    else
                        total_pos = (obj.current_max_index + 1)^N - 1;
                        T = zeros(dimension_vector(d, 2 * (N + 1), [1, 1]));
                    end


                    %generate all combinations of internam indices
                    parfor i = 0:total_pos

                        if cyclic == 1
                            full_vect = [encode_index_array(i, N + 1, obj.current_max_index + 1); 0];
                            full_vect(end) = full_vect(1);
                        else
                            full_vect = [virtual_level; encode_index_array(i, N, obj.current_max_index + 1); virtual_level];
                        end


                        correct_index_set = 1;

                        O_tensors = cell(1, N+1);
                        %Take correct O's for contraction with ncon
                        for j = 1:N + 1
                            O_j = O{full_vect(j)+1, full_vect(j + 1)+1};
                            if length(O_j) == 0
                                correct_index_set = 0;
                                break;
                            end
                            O_tensors{j} = O_j;
                        end


                        if correct_index_set == 1
                            legs = cell(1, N+1);


                            for t = 1:N
                                legs{t} = zeros(1, 4);
                            end

                            legs{1}(1) = -1; %first index
                            legs{N+1}(4) = -(2 * (N + 1) + 2); % last one

                            for t = 1:N %assign all indices to sum over
                                legs{t}(4) = t;
                                legs{t+1}(1) = t;
                            end

                            for t = 1:N + 1 %assign final places for i_n and j_n
                                legs{t}(2) = -(t + 1) + cyclic;
                                legs{t}(3) = -(N + 1 + t + 1) + cyclic;
                            end

                            if cyclic == 1
                                legs{1}(1) = N + 1;
                                legs{end}(4) = N + 1;
                            end

                            T_j = ncon(O_tensors, legs);
                            T = T + T_j;
                        end
                    end

                case "matrix"

                    MPO = obj.MPO_cell;


                    if cyclic == 1
                        M = N + 1;
                        tensors = cell(1, M);
                        tensors(1:end) = {MPO};

                        leg_list = cell(1, M);

                        for i = 1:M
                            leg_list{i} = [i - 1, -(i), -(i + M), i];
                        end

                        leg_list{1}(1) = leg_list{end}(4);

                        T = ncon(tensors, leg_list);

                    else

                        M = N + 3;

                        tensors = cell(1, M);
                        tensors(2:end-1) = {MPO};
                        tensors{1} = obj.left;
                        tensors{M} = obj.right;


                        leg_list = cell(1, M);
                        leg_list{1} = [-1, 1];
                        leg_list{M} = [M - 1, -2 * (M - 1)];

                        for i = 2:M - 1
                            leg_list{i} = [i - 1, -(i), -(i + M - 2), i];
                        end

                        T = ncon(tensors, leg_list);

                    end

                otherwise
                    error("unknown type")
            end


            if matrix == 1
                T = reshape(T, [obj.dim^(N + 1), obj.dim^(N + 1)]);
                return
            end


        end

        function [left_i, left_dim] = get_L(obj, M, U_cell)

            cellsource = 0;

            if nargin < 3
                cellsource = 1;
            end

            d = obj.dim;

            if M ~= 0

                left_list = cell(1, M);
                contract_list = cell(1, M);

                for i = 1:M

                    if cellsource == 1
                        left_list{i} = obj.MPO_cell{i, i+1};
                    else
                        left_list{i} = U_cell{i};
                    end


                    contract_list{i} = [i, -(2 * i), -(2 * i + 1), i + 1];

                end

                contract_list{1}(1) = -1;
                contract_list{end}(4) = -(2 * M + 2);


                %auto detect size
                left_i = reshape(ncon(left_list, contract_list), d^(2 * M), []);

                left_dim = size(left_i, 2);

            else
                left_dim = 1;
                left_i = [1];
            end

        end

        function [left_inv, left_dim] = get_L_inv_MPO(obj, M)
            if M ~= 0

                left_list = cell(1, M);
                contract_list = cell(1, M);

                for i = 1:M
                    left_list{i} = obj.inverse_MPO_left{i};
                    contract_list{i} = [i, i - 1, -(2 * i), -(2 * i + 1)];
                end


                contract_list{1}(2) = -(2 * M + 2);
                contract_list{end}(1) = -1;


                left_inv = reshape(ncon(left_list, contract_list), [], obj.dim^(2 * M));
                left_dim = size(left_inv, 1);

            else
                left_inv = [1];
                left_dim = 1;
            end

        end

        function [right_inv, right_dim] = get_R_inv_MPO(obj, M)
            if M ~= 0

                right_list = cell(1, M);
                contract_list = cell(1, M);

                for i = 1:M
                    j = M - i + 1;
                    right_list{i} = obj.inverse_MPO_right{i};
                    contract_list{i} = [-(2 * j), -(2 * j + 1), i - 1, i];

                end

                contract_list{1}(3) = -1;
                contract_list{end}(4) = -(2 * M + 2);

                right_inv = reshape(ncon(right_list, contract_list), obj.dim^(2 * M), []);


                right_dim = size(right_inv, 2);

            else
                right_inv = [1];
                right_dim = 1;
            end

        end

        function [right_i, right_dim] = get_R(obj, M, V_cell)

            cellsource = 0;

            if nargin < 3
                cellsource = 1;
            end

            d = obj.dim;

            if M ~= 0


                right_list = cell(1, M);
                contract_list = cell(1, M);

                for i = 1:M

                    if cellsource == 1
                        right_list{end-i+1} = obj.MPO_cell{i+1, i};
                    else
                        right_list{end-i+1} = V_cell{i};
                    end


                    contract_list{i} = [i, -(2 * i), -(2 * i + 1), i + 1];

                end

                contract_list{1}(1) = -1;
                contract_list{end}(4) = -(2 * M + 2);


                right_i = reshape(ncon(right_list, contract_list), [], d^(2 * M));

                right_dim = size(right_i, 1);

            else
                right_dim = 1;
                right_i = [1];
            end

        end

        function [y, left_dim] = invert_left(obj, x, M, U_cell)

            if obj.inverse_MPO
                %recollect inverses from memory
                [left_inv, left_dim] = get_L_inv_MPO(obj, M);
                y = ncon({left_inv, x}, {[-1, 1], [1, -2]});
            else
                if nargin < 4
                    [left_i, left_dim] = get_L(obj, M);
                else
                    [left_i, left_dim] = get_L(obj, M, U_cell);
                end
                y = lsqminnorm(left_i, x);
            end

        end

        function [y, right_dim] = invert_right(obj, x, M, V_cell)
            %
            if obj.inverse_MPO
                %recollect inverses from memory
                [right_inv, right_dim] = get_R_inv_MPO(obj, M);
                y = ncon({x, right_inv}, {[-1, 1], [1, -2]});
            else
                if nargin < 4
                    [right_i, right_dim] = get_R(obj, M);
                else
                    [right_i, right_dim] = get_R(obj, M, V_cell);
                end
                y = lsqminnorm(right_i', x')';
            end
        end

        function obj = register_left_inv_mpo(obj, N, U_inv, sigma)

            if obj.inverse_MPO

                ldim = size(U_inv, 2) / (obj.dim^2);

                if nargin < 4
                    rdim = size(U_inv, 1);
                    inv = U_inv;
                else
                    rdim = size(sigma, 2);
                    invSigma = diag(sigma).^-1;
                    inv = diag(invSigma) * U_inv;
                end


                invMPO = reshape(inv, [rdim, ldim, obj.dim, obj.dim]);

                obj.inverse_MPO_left{N+1, 1} = invMPO;
            end
        end

        function obj = register_right_inv_mpo(obj, N, V_inv, sigma)
            if obj.inverse_MPO

                rdim = size(V_inv, 1) / obj.dim^2;

                if nargin < 4
                    ldim = size(V_inv, 2);
                    inv = V_inv;
                else
                    ldim = size(sigma, 1);
                    invSigma = diag(sigma).^-1;
                    inv = V_inv * diag(invSigma);
                end

        
                invMPO = reshape(inv, [obj.dim, obj.dim, rdim, ldim]);

                obj.inverse_MPO_right{N+1, 1} = invMPO;
            end
        end

        function [O, normalisation_factor] = change_normalisation(O, normalisation_factor, change)
            max_dim = size(O, 1);

            for i = 1:max_dim

                for j = 1:max_dim
                    O{i, j} = O{i, j} * change;
                end
            end

            normalisation_factor = normalisation_factor / change;

        end
    end
end

%helper functions


function p = dimension_vector(d, n, leftright)
%helper function to create a 1xn vector  [ left,d,d,..,d,right]
%if left/right are not supplied/0, this is omitted

if nargin < 3
    p = zeros(1, n);
    p = p + d;
    return

else

    p = zeros(1, n+2);
    p = p + d;
    p(1) = leftright(1);
    p(end) = leftright(2);
end

end

function y = encode_index_array(n, len, d)
% this takes a single number and splits it into the composite indices. The
i = 1;
y = zeros(len, 1);

while n ~= 0
    [n, r] = readOne(n, d);
    y(i) = r;
    i = i + 1;
end

end

function [s, r] = readOne(s, d)
r = rem(s, d);
s = (s - r) / d;
end

function p = site_ordering_permute(n)
% changes from |left i1 i2 ... j1 j2.. right> to |left i1 j1 i2 j2 ...
% right>
p = zeros(2*n+2, 1);
p(1) = 1;
p(2*n+2) = 2 * n + 2;
for i = 1:n
    p(2*i) = i + 1;
end
for i = 1:n
    p(2*i+1) = n + i + 1;
end
end

%just the element wise 2 norm
function norm = tensor_norm(X)
v = reshape(X, [], 1);
%N=length(v);
norm = sqrt(sum(v.^2));
end

function z = geomSum(x, n)
z = (x^(n + 1) - 1) / (x - 1);
end
