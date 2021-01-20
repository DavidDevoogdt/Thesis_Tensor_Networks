%for tensors:       ( beta)
%            (alpha)-- O -- (gamma)  = O(i,j,alpha,beta,gamma,delta)
%                    (delta)             1 2   3     4    5     6
% for tensors containig multiple ij's: O numbered from left to right and
% for a given vertical pos from up till down

classdef PEPO

    properties
        dim
        H_1_tensor
        H_2_tensor
        PEPO_cell
        type
        nf%normalisation factor
        max_index
        testing
        visualise
        virtual_level_sizes_horiz
        virtual_level_sizes_vert
        PEPO_matrix
        current_max_index
        numopts
        cycleopts
        cycle_index
        boundary_matrix_x
        boundary_matrix_y
        order
    end

    methods

        function obj = PEPO(d, H_1_tensor, H_2_tensor, order, type, opts)
            numopts.numbered = 1;
            obj.numopts = numopts;

            cycleopts.numbered = 1;
            cycleopts.v_cyclic = 0;
            cycleopts.h_cyclic = 1;
            obj.cycleopts = cycleopts;

            obj.dim = d;
            obj.H_1_tensor = H_1_tensor;
            obj.H_2_tensor = H_2_tensor;

            if mod(order, 2)
                max_index = (order - 1) / 2;
            else
                max_index = order / 2;
            end

            obj.order = order;

            obj.PEPO_cell = cell(max_index + 1, max_index + 1, max_index + 1, max_index + 1);
            obj.boundary_matrix_x = cell(max_index + 1, max_index + 1);
            obj.boundary_matrix_y = cell(max_index + 1, max_index + 1);
            obj.boundary_matrix_x{1, 1} = reshape([1], 1, 1);
            obj.boundary_matrix_y{1, 1} = reshape([1], 1, 1);

            obj.type = type;
            obj.max_index = max_index;

            obj.cycle_index = Inf;

            p = inputParser;
            addParameter(p, 'testing', 0)
            addParameter(p, 'visualise', 0)
            parse(p, opts)

            obj.testing = p.Results.testing;
            obj.visualise = p.Results.visualise;

            %calculate ln of normailisation fac

            map = create_map(1:2, numopts);
            [~, nf2] = H_exp(obj, map, 0, true);

            obj.nf = nf2;

            obj = obj.makePEPO();
        end

        function obj = makePEPO(obj)
            d = obj.dim;

            %%%%%%%%%%single site

            obj.PEPO_cell{1, 1, 1, 1} = reshape(eye(d) / exp(obj.nf), [d, d, 1, 1, 1, 1]);

            %             [map,~] = create_map(1:2,obj.numopts);
            %             [target0,ln_prefact_out] = obj.H_exp(map,obj.nf,false  );
            %
            %             target0_site = reshape( permute(target0, site_ordering_permute(map.N) ), dimension_vector(d^2,map.N)  );
            %
            %             pattern0 = {[0,0,0,0]};
            %
            %             con_cells0 = get_valid_contractions(obj,map, struct('max_index', 0,'pattern',{pattern0} ) );
            %
            %             x_cell0 = obj.solve_non_lin( pattern0,{map},{target0_site},{con_cells0},struct('Gradient',false), ln_prefact_out );
            %
            %             O_0000 = x_cell0{1};
            %             %O_0000 = (x_cell0{1}+ x_cell0{1}')/2 ;   %expm( 0*obj.H_1_tensor ); % eye(d);%expm( 0*obj.H_1_tensor );
            %             %obj.nf = trace(O_0000);
            %
            %             obj.nf = -log( svds( O_0000,1)  )+obj.nf;
            %
            %             diff = O_0000-obj.PEPO_cell{1,1,1,1};
            %            obj.PEPO_cell{1,1,1,1} = reshape(  O_0000/exp(obj.nf) , [d,d,1,1,1,1] ) ;
            %%%%%%%%%% identity

            function [map, con_cells_b, target, res_target, ln_prefact_out, rank_x] = double(n, strict, non_lin, ln_prefact)

                obj.current_max_index = n / 2;

                [map, ~] = create_map(1:n, obj.numopts);
                [target, ln_prefact_out] = H_exp(obj, map, ln_prefact, true);
                target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(d^2, map.N));

                target = reshape(target, [d^map.N, d^map.N]);

                m = n / 2;

                pattern = {[m - 1, 0, m, 0], [m, 0, m - 1, 0]};

                con_cells = get_valid_contractions(obj, map, struct('max_index', obj.current_max_index, 'pattern', {pattern}));
                [con_cells_a, con_cells_b] = split_con_cells(map, con_cells);

                mul_factor = exp(ln_prefact_out - obj.nf);

                cc = con_cells;

                [x_cell, res_target, rank_x, mask] = solve_lin(obj, pattern, map, cc, target_site, ln_prefact_out);
                %
                %                 if sum(mask) ~= 0
                %
                %                     err_before = obj.calculate_error(1:n+1,obj.numopts) ;
                %
                %                     obj.PEPO_cell{m,1,m+1,1} = x_cell{1}*mul_factor;%right
                %                     obj.PEPO_cell{m+1,1,m,1} = x_cell{2}*mul_factor;%left
                %
                %                     err_after= obj.calculate_error(1:n+1,obj.numopts) ;
                %
                %
                %                     if err_after>err_before
                %                         mul_factor = 0;
                %                         rank_x = 0;
                %                     end
                %                 end
                %

                if rank_x ~= 0

                    obj.PEPO_cell{m, 1, m + 1, 1} = x_cell{1} * mul_factor; %right
                    obj.PEPO_cell{m + 1, 1, m, 1} = x_cell{2} * mul_factor; %left
                end

                %
                res_target = ipermute(reshape(res_target, dimension_vector(d, 2 * map.N)), site_ordering_permute(map.N));

                target = reshape(target, [d^map.N, d^map.N]);
                res_target = reshape(res_target, [d^map.N, d^map.N]);

                if obj.testing == 1
                    err = obj.calculate_error(1:n, obj.numopts)
                end

            end

            function [map, con_cells_b, target, res_target, ln_prefact_out, rank_x] = single(n, strict, non_lin, ln_prefact)
                %%%%%%%%%%%%%%%%  determine 1--|--1
                [map, ~] = create_map(1:n, obj.numopts);

                m = (n - 1) / 2;

                %d2 = d^4;

                [target, ln_prefact_out] = H_exp(obj, map, ln_prefact, true);

                target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(d^2, map.N));

                target = reshape(target, [d^map.N, d^map.N]);

                pattern = {[m, 0, m, 0]};

                con_cells = get_valid_contractions(obj, map, struct('max_index', obj.current_max_index, 'pattern', {pattern}));
                [con_cells_a, con_cells_b] = split_con_cells(map, con_cells);

                cc = con_cells;

                [x_cell, res_target, rank_x, mask] = solve_lin(obj, pattern, map, cc, target_site, ln_prefact_out);

                mul_factor = exp(ln_prefact_out - obj.nf);

                %
                %                  if sum(mask) ~= 0
                %
                %                     err_before = obj.calculate_error(1:n+1,obj.numopts) ;
                %
                %                     obj.PEPO_cell{m+1,1,m+1,1} =  x_cell{1}*mul_factor;
                %
                %                     err_after= obj.calculate_error(1:n+1,obj.numopts) ;
                %
                %
                %                     if err_after>err_before
                %                         mul_factor = 0;
                %                         rank_x = 0;
                %                     end
                %                   end

                obj.PEPO_cell{m + 1, 1, m + 1, 1} = x_cell{1} * mul_factor;
                %
                %                 [map2,~] = create_map(1:n+1,obj.cycleopts);
                %
                %                 [target2,~] = obj.H_exp(map2,ln_prefact_out,false  );
                %
                %
                %                 target2_site = reshape( permute(target2, site_ordering_permute(map2.N) ), dimension_vector(d^2,map2.N)  );
                %
                %                 con_cells2 = get_valid_contractions(obj,map2, struct('max_index', obj.current_max_index,'pattern',{pattern}  ));
                %
                %                 [con_cells2, target2_site] = obj.optimize_con_cells({map2}, {con_cells2} , pattern,{target2_site}, ln_prefact_out );
                %
                %
                %                 con_cells = { cc, con_cells2{1} };
                %                 targets = {target_site, target2_site{1}};
                %                 maps = {map, map2};
                %
                %
                %                 x_cell_2 = obj.solve_non_lin(pattern,maps,targets,con_cells, struct(),ln_prefact_out  );
                %
                %
                %
                %
                %
                %                 obj.PEPO_cell{m+1,1,m+1,1} =  x_cell_2{1}*mul_factor;
                %

                res_target = ipermute(reshape(res_target, dimension_vector(d, 2 * map.N)), site_ordering_permute(map.N));

                res_target = reshape(res_target, [d^map.N, d^map.N]);

                if obj.testing == 1
                    err = obj.calculate_error(1:n, obj.numopts)
                end

            end

            %%%%%%%%%%%%%% 0--|--1--|--0 and all other veriants
            obj.current_max_index = 1;
            obj.virtual_level_sizes_horiz = [1];
            obj.virtual_level_sizes_vert = [1];

            ln_prefact = obj.nf;

            for n = 2:obj.order
                strict = false;
                non_lin = false;

                if mod(n, 2) == 1

                    [~, con_cells_b, target, res_target, ln_prefact, rank_x] = single(n, strict, non_lin, ln_prefact);

                    e1 = obj.calculate_error(1:n + 1, obj.numopts);

                    %                 e1 = svds(res_target,1);
                    e2 = svds(target, 1);
                    %                 rel_residual = e1/e2;
                    fprintf("n=%d residual = %.4e  tar %.4e d_nf %.4e  \n", n, e1, e2, exp(ln_prefact - obj.nf));
                    %
                    if rank_x == 0
                        obj.max_index = obj.current_max_index;

                        break;
                    end

                else
                    obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^(n)];
                    obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^(n)];

                    obj.current_max_index = n / 2;

                    [~, con_cells_b, target, res_target, ln_prefact, rank_x] = double(n, strict, non_lin, ln_prefact);

                    e1 = obj.calculate_error(1:n + 1, obj.numopts);

                    %                 e1 = svds(res_target,1);
                    e2 = svds(target, 1);
                    %                 rel_residual = e1/e2;
                    fprintf("n=%d residual = %.4e, tar %.4e d_nf %.4e \n", n, e1, e2, exp(ln_prefact - obj.nf))

                    if rank_x == 0

                        obj.virtual_level_sizes_horiz = obj.virtual_level_sizes_horiz(1:end - 1);
                        obj.virtual_level_sizes_vert = obj.virtual_level_sizes_vert(1:end - 1);

                        obj.current_max_index = obj.current_max_index - 1;
                        obj.max_index = obj.current_max_index;

                        break;
                    end

                end

            end

            %%%%%%%%%%%

            obj = obj.cell2matrix(); %save matrix form

            %if obj.testing == 1
            %                  max(reshape(abs(obj.PEPO_cell{1,1,1,1}),[],1))
            %                  max(reshape(abs(obj.PEPO_cell{2,1,2,1}),[],1))
            %                  max(reshape(abs(obj.PEPO_cell{1,1,2,1}),[],1))
            %                  max(reshape(abs(obj.PEPO_cell{2,1,3,1}),[],1))
            %                  max(reshape(abs(obj.PEPO_cell{3,1,3,1}),[],1))
            %                  obj.calculate_error(1:6,obj.cycleopts)
            %                  obj.calculate_error(1:7,obj.cycleopts)
            %                  obj.calculate_error(1:8,obj.cycleopts)
            %
            %                  obj.calculate_error(1:6,obj.numopts)
            %                  obj.calculate_error(1:7,obj.numopts)
            %                  obj.calculate_error(1:8,obj.numopts)

            %end

            hsize = sum(obj.virtual_level_sizes_horiz);
            T = reshape(obj.PEPO_matrix(1, 1, :, 1, :, 1), [hsize, hsize]) * exp(obj.nf);

            T2 = T;
            left = zeros(hsize, 1);
            left(1) = 1;

            S = svds(T2) / exp(ln_prefact);
        end

        function obj = makePEPO3(obj)
            d = obj.dim;

            %todo do this in code
            obj.virtual_level_sizes_horiz = [d^0, d^2, d^2];
            obj.virtual_level_sizes_vert = [d^0, d^2, d^2];

            %%%%%%%%%%single site
            O_0000 = eye(d); %expm( 0*obj.H_1_tensor );
            %obj.nf = trace(O_0000);

            obj.PEPO_cell{1, 1, 1, 1} = reshape(O_0000 / exp(obj.nf), [d, d, 1, 1, 1, 1]);

            %%%%%%%%%%%%%% 0--|--1--|--0 and all other veriants
            obj.current_max_index = 0;

            part = obj.get_middle_part(...
                {[], [], [], []}, [1, 2], 0);

            [U, S, V] = svd(reshape(part, d^2, d^2));

            sqrt_S = diag(diag(S).^0.5);

            block_l = permute(reshape(U * sqrt_S, [1, d, d, d^2]), [2, 3, 1, 4]);
            block_r = permute(reshape(sqrt_S * V', [d^2, d, d, 1]), [2, 3, 1, 4]);

            obj.PEPO_cell{1, 1, 2, 1} = reshape(block_l, [d, d, 1, 1, d^2, 1]); %right
            obj.PEPO_cell{1, 1, 1, 2} = reshape(block_l, [d, d, 1, 1, 1, d^2]); %down

            obj.PEPO_cell{2, 1, 1, 1} = reshape(block_r, [d, d, d^2, 1, 1, 1]); %left
            obj.PEPO_cell{1, 2, 1, 1} = reshape(block_r, [d, d, 1, d^2, 1, 1]); %up

            %%%%%%%%%%%%%%%%%create 0--|--1--|--1--|--0 and variants
            if obj.max_index >= 1

                obj.current_max_index = 1;

                obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz, d^4]
                obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert, d^4]

                %solve with new solver;
                map = create_map([1, 2, 3], obj.numopts);

                Tensor = obj.H_exp(map, obj.nf) - ...
                    obj.contract_network(map, struct('max_index', obj.current_max_index));

                Tensor_site = reshape(permute(Tensor, site_ordering_permute(map.N)), ...
                    [d^2, d^2, d^2]);

                opts.tol = 1e-13;
                opts.maxit = 1;
                opts.print_level = 1;
                opts.optim = [2]; %only optimize 2;
                opts.get_elem_num = [1, 2, 3];
                opts.solve_type = {'', 'matrix_inv', ''};

                elem_list = cell(3, 1);
                elem_list{1} = obj.PEPO_cell{1, 1, 2, 1};
                elem_list{2} = rand(d, d, d^2, 1, d^2, 1);
                elem_list{3} = obj.PEPO_cell{2, 1, 1, 1};

                [elems, err, A] = tensor_ring(elem_list, map, Tensor_site, opts);

                %
                block_11 = elems{2};
                %block_11_bis = obj.get_middle_part(
                %{[1,2],[],[],[2;3]},[1,2;0,3] ); same

                %6E-8
                %             obj.calculate_error( create_map([1 2 3],1)) %  0
                %             obj.calculate_error( create_map([1 2;0 3],1)) % 0
                %             obj.calculate_error( create_map([1 0; 2 3],1)) %2E-16
                %             obj.calculate_error( create_map([1; 2 ;3;],1)) %  0

                %al same block because up and left blocks are the same
                obj.PEPO_cell{2, 1, 2, 1} = reshape(block_11, [d, d, d^2, 1, d^2, 1]);
                obj.PEPO_cell{2, 1, 1, 2} = reshape(block_11, [d, d, d^2, 1, 1, d^2]);
                obj.PEPO_cell{1, 2, 2, 1} = reshape(block_11, [d, d, 1, d^2, d^2, 1]);
                obj.PEPO_cell{1, 2, 1, 2} = reshape(block_11, [d, d, 1, d^2, 1, d^2]);
                %
                %

                if obj.testing == 1
                    obj.calculate_error(create_map([1 2 3], obj.numopts))
                    obj.calculate_error(create_map([1 2; 0 3], obj.numopts))
                    obj.calculate_error(create_map([1 0; 2 3], obj.numopts))
                    obj.calculate_error(create_map([1; 2; 3; ], obj.numopts))
                end

                %make cyclic properties better with null space vectors
                opts_h_cyclic_numbered.numbered = 1;
                opts_h_cyclic_numbered.v_cyclic = 0;
                opts_h_cyclic_numbered.h_cyclic = 1;

                map_c = create_map([1, 2, 3], opts_h_cyclic_numbered); %cyclic ring

                Tensor_c = obj.H_exp(map_c, obj.nf) - ...
                    obj.contract_network(map_c, struct('max_index', obj.current_max_index));

                Tensor_site_c = reshape(permute(Tensor_c, site_ordering_permute(map.N)), ...
                    [d^2, d^2, d^2]);

                opts.tol = 1e-13;
                opts.maxit = 1;
                opts.print_level = 1;
                opts.get_elem_num = [1; 1; 1];
                opts.solve_type = {'fsolve', '', ''};
                opts.null_space

                elem_list{1} = rand(d, d, d^2, 1, d^2, 1);
                [elems, err, A] = tensor_ring(elem_list, map, Tensor_site, opts);

                obj.PEPO_cell{3, 1, 3, 1} = elems{1};

                if obj.testing == 1
                    obj.calculate_error(create_map([1 2 3], obj.numopts))%no improvement
                    obj.calculate_error(create_map([1 2 3], opts_h_cyclic_numbered))
                    obj.calculate_error(create_map([1 2 3 4 5 6], opts_h_cyclic_numbered))
                end

            end

            %%%%%%%%%%%

            obj = obj.cell2matrix(); %save matrix form

        end

        function [A, B, G1, lambda1] = vumps(obj, chimax)

            %todo check these params
            opts.charges = 'regular';
            opts.dynamical = 'off';
            opts.dyncharges = 0;
            opts.schmidtcut = 1e-10;
            opts.chimax = 350;
            %opts.disp='iter';
            opts.disp = 'none';
            opts.tolmax = 1e-4; %1e-4
            opts.tolfactor = 1e4;
            opts.minit = 1;
            opts.dyniter = 5;
            opts.truncate = 0;
            opts.method = 'vumps';
            opts.save = 0;

            %opts.method = 'qr';

            opts.plot = 'on';
            opts.maxit = 1000;
            opts.tolfixed = 1e-12;

            %put into vumps format

            T = obj.PEPO_matrix;

            hdim = size(T, 3);
            vdim = size(T, 4);

            %upper vumps zipper
            M = ncon({T}, {[1, 1, -1, -2, -3, -4]});
            %lower vumps zipper

            o.legs = 4;
            o.group = 'none';
            o.dims = size(M);
            o.var = M;

            O.type = 'mpo';
            O.mpo = o;

            [A, G1, lambda1, ~, ~] = Vumps(O, chimax, [], opts);

            %correct estimate for inversion sym?

            GL = G1{1}; GR = G1{2}; Ac = A{4};

            m.legs = 4;
            m.group = 'none';
            m.dims = size(M);
            m.var = M;

            opts.krylovdim = 100; opts.tol = 1e-14;
            %opts.disp='iter-detailed'; %opts.reorth='force';

            function x = vumps_under(x)
                x = TensorContract({x, GL, m, GR}, {[1, 2, 5], [1, 3, -1], [-2, 4, 2, 3], [-3, 4, 5]});
            end

            [B, lambda2, err2] = TensorEigs(@(x) vumps_under(x), TensorConj(Ac), 1, 'lm', opts);

            %[B_l,C_l,~]=TensorDecRight(Bc,'polar');

            %[B_r,C_r,~]=TensorDecLeft(Bc,'polar');

            %             accopts.method = 'qr';
            %             accopts.tol = 1e-16;
            %
            %             [B,err]=VumpsSolveACC(B2,B_c,accopts);

        end

        function [mag, corr_length, delta] = get_expectation (obj, X, chimax)
            [A, B, G, lambda] = vumps(obj, chimax);

            T = obj.PEPO_matrix;

            M = ncon({T}, {[1, 1, -1, -2, -3, -4]});

            m.legs = 4;
            m.group = 'none';
            m.dims = size(M);
            m.var = M;

            O = ncon ({T, X}, {[1, 2, -1, -2, -3, -4], [1, 2]});
            o.legs = 4;
            o.group = 'none';
            o.dims = size(O);
            o.var = O;

            %transfereigs edited !

            GL = G{1}; GR = G{2}; Ac = A{4};

            %should be lambda?
            [x, ~] = TensorContract({B, GL, Ac, m, GR}, ...
                {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3], [8, 7, 6]});

            [y, ~] = TensorContract({B, GL, Ac, o, GR}, ...
                {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3], [8, 7, 6]});

            mag = y / x;

            function x = transfer_up(x)
                [x, ~] = TensorContract({GL, x, m, GR}, ...
                    {[-1, 3, 4], [4, 5, 8], [5, 7, -2, 3], [8, 7, -3]});
            end

            opts.krylovdim = 100; opts.tol = 1e-14;

            [rho, f] = TensorEigs(@(x) transfer_up(x), A{4}, 5, 'lm', opts);

            eps_i = -log(abs(f));
            corr_length = eps_i(1 + 1) - eps_i(1);

            delta = eps_i(4 + 1) - eps_i(2 + 1); %same as https://arxiv.org/pdf/1907.08603.pdf
        end

    end

end
