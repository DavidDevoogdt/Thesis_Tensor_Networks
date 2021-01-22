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

            e1 = calculate_error(obj, 1:n + 1, obj.numopts);

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

            e1 = calculate_error(obj, 1:n + 1, obj.numopts);

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

    obj = cell2matrix(obj); %save matrix form

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
