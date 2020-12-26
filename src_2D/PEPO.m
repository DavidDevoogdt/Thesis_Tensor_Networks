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
        nf %normalisation factor
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
    end
    
    
        
    methods
        function obj = PEPO(d,H_1_tensor,H_2_tensor,max_index,type,opts)
            numopts.numbered = 1;
            obj.numopts = numopts;
            
             cycleopts.numbered=1;
             cycleopts.v_cyclic=0;
             cycleopts.h_cyclic=1;
            obj.cycleopts = cycleopts;
            
            obj.dim = d;
            obj.H_1_tensor = H_1_tensor;
            obj.H_2_tensor = H_2_tensor;
            obj.PEPO_cell = cell( max_index+1,max_index+1,max_index+1,max_index+1 );
            obj.boundary_matrix_x = cell( max_index+1,max_index+1);
            obj.boundary_matrix_y = cell( max_index+1,max_index+1);
            obj.boundary_matrix_x{1,1} = reshape([1],1,1);
            obj.boundary_matrix_y{1,1} = reshape([1],1,1);
            
            obj.type = type;
            obj.max_index = max_index;
            
            obj.cycle_index = Inf;
            
            p = inputParser;
            addParameter(p,'testing',0)
            addParameter(p,'visualise',0)
            parse(p,opts)
            
            obj.testing= p.Results.testing;
            obj.visualise = p.Results.visualise;
            
               %nf calculated for square
            map = obj.create_map( [1,1;1,1]);

            Hexp = reshape(  obj.H_exp(map) , d^4,d^4);
            
            S = sum(abs(svds(Hexp)));
            obj.nf = S.^(1/4);
            
            
            obj = obj.makePEPO();
        end
        
        function obj = makePEPO3(obj)
            d = obj.dim;
            
            %todo do this in code
            obj.virtual_level_sizes_horiz = [d^0,d^2,d^2];
            obj.virtual_level_sizes_vert = [d^0,d^2,d^2];

            
            %%%%%%%%%%single site
            O_0000 = eye(d);%expm( 0*obj.H_1_tensor );
            %obj.nf = trace(O_0000);
            
            obj.PEPO_cell{1,1,1,1} = reshape(  O_0000/obj.nf , [d,d,1,1,1,1] ) ;

            %%%%%%%%%%%%%% 0--|--1--|--0 and all other veriants
            obj.current_max_index = 0;
            
            
            part = obj.get_middle_part(...
                {[],[],[],[]},[1,2],0);


            [U,S,V] = svd(  reshape(part,d^2,d^2) );
            
            sqrt_S = diag(diag(S).^0.5);
            
            block_l = permute( reshape(U*sqrt_S, [1,d,d,d^2]), [2,3,1,4]);
            block_r = permute( reshape(sqrt_S*V', [d^2,d,d,1]), [2,3,1,4]); 
            

            obj.PEPO_cell{1,1,2,1} =reshape(block_l, [d,d,1,1,d^2,1]);%right
            obj.PEPO_cell{1,1,1,2} =reshape(block_l, [d,d,1,1,1,d^2]);%down

            obj.PEPO_cell{2,1,1,1} =reshape(block_r, [d,d,d^2,1,1,1]);%left
            obj.PEPO_cell{1,2,1,1} =reshape(block_r, [d,d,1,d^2,1,1]);%up
            

            
            %%%%%%%%%%%%%%%%%create 0--|--1--|--1--|--0 and variants
           if obj.max_index >=1 
            
               obj.current_max_index=1;
               
               
                obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz,d^4] 
                obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert,d^4] 
               
                %solve with new solver;
                map = PEPO.create_map([1,2,3],obj.numopts); 
              
                Tensor = obj.H_exp(map,obj.nf)-...
                    obj.contract_network(map,struct('max_index', obj.current_max_index));

                Tensor_site = reshape(  permute(Tensor, site_ordering_permute(map.N)),...
                            [d^2,d^2,d^2] );        
                

                opts.tol = 1e-13;
                opts.maxit = 1;
                opts.print_level = 1;
                opts.optim = [2]; %only optimize 2;
                opts.get_elem_num = [1,2,3];
                opts.solve_type = {'','matrix_inv',''};
                
                elem_list = cell(3,1  ) ;
                elem_list{1} =   obj.PEPO_cell{1,1,2,1};
                elem_list{2} =   rand( d,d, d^2,1,d^2,1);
                elem_list{3} =   obj.PEPO_cell{2,1,1,1};
                
                [elems,err,A] = tensor_ring(elem_list,map,Tensor_site, opts);
                
%                 
                block_11 = elems{2};
                %block_11_bis = obj.get_middle_part(
                %{[1,2],[],[],[2;3]},[1,2;0,3] ); same 

                %6E-8
    %             obj.calculate_error( PEPO.create_map([1 2 3],1)) %  0
    %             obj.calculate_error( PEPO.create_map([1 2;0 3],1)) % 0
    %             obj.calculate_error( PEPO.create_map([1 0; 2 3],1)) %2E-16
    %             obj.calculate_error( PEPO.create_map([1; 2 ;3;],1)) %  0

                %al same block because up and left blocks are the same
                obj.PEPO_cell{2,1,2,1}= reshape(block_11,[d,d,d^2,1,d^2,1]); 
                obj.PEPO_cell{2,1,1,2}= reshape(block_11,[d,d,d^2,1,1,d^2]); 
                obj.PEPO_cell{1,2,2,1}= reshape(block_11,[d,d,1,d^2,d^2,1]); 
                obj.PEPO_cell{1,2,1,2}= reshape(block_11,[d,d,1,d^2,1,d^2]); 
% 
%                 

                if obj.testing ==1
                    obj.calculate_error( PEPO.create_map([1 2 3],obj.numopts)) 
                    obj.calculate_error( PEPO.create_map([1 2;0 3],obj.numopts)) 
                    obj.calculate_error( PEPO.create_map([1 0; 2 3],obj.numopts)) 
                    obj.calculate_error( PEPO.create_map([1; 2 ;3;],obj.numopts)) 
                end
               
               
               %make cyclic properties better with null space vectors
               opts_h_cyclic_numbered.numbered=1;
               opts_h_cyclic_numbered.v_cyclic=0;
               opts_h_cyclic_numbered.h_cyclic=1;
                
               map_c = PEPO.create_map([1,2,3],opts_h_cyclic_numbered); %cyclic ring
              
                Tensor_c = obj.H_exp(map_c,obj.nf)-...
                    obj.contract_network(map_c,struct('max_index', obj.current_max_index));

                Tensor_site_c = reshape(  permute(Tensor_c, site_ordering_permute(map.N)),...
                            [d^2,d^2,d^2] );        

                opts.tol = 1e-13;
                opts.maxit = 1;
                opts.print_level = 1;
                opts.get_elem_num = [1;1;1];
                opts.solve_type = {'fsolve','',''};
                opts.null_space
                
               
                elem_list{1} = rand(  d,d, d^2,1,d^2,1);
                [elems,err,A] = tensor_ring(elem_list,map,Tensor_site, opts);
          
               
                obj.PEPO_cell{3,1,3,1} = elems{1} ;
                
                if obj.testing ==1
                    obj.calculate_error( PEPO.create_map([1 2 3],obj.numopts)) %no improvement
                    obj.calculate_error( PEPO.create_map([1 2 3],opts_h_cyclic_numbered))
                    obj.calculate_error( PEPO.create_map([1 2 3 4 5 6],opts_h_cyclic_numbered))
                end
                
           end
            
            
            
            %%%%%%%%%%%

            obj = obj.cell2matrix() ; %save matrix form
           
        end
%       
        function obj = makePEPO(obj)
            d = obj.dim;
            
            %todo do this in code
            obj.virtual_level_sizes_horiz = [d^0,d^2];
            obj.virtual_level_sizes_vert = [d^0,d^2];

            
            %%%%%%%%%%single site
            O_0000 =expm( 0*obj.H_1_tensor ); % eye(d);%expm( 0*obj.H_1_tensor );
            %obj.nf = trace(O_0000);
            
            obj.PEPO_cell{1,1,1,1} = reshape(  O_0000/obj.nf , [d,d,1,1,1,1] ) ;

            
            
            %%%%%%%%%%%%%% 0--|--1--|--0 and all other veriants
            obj.current_max_index = 0;
            
            
            part = obj.get_middle_part(...
                {[],[],[],[]},[1,2],0);


            [U,S,V] = svd(  reshape(part,d^2,d^2) );
            
            n_extra = 4;
            d1 = d^2+n_extra;
            
            [U,S,V] = expand_svd(U,S,V,n_extra);
            
            sqrt_S = diag(diag(S).^0.5);
            
            block_l = permute( reshape(U*sqrt_S, [1,d,d,d1]), [2,3,1,4]);
            block_r = permute( reshape(sqrt_S*V', [d1,d,d,1]), [2,3,1,4]); 
            

            obj.PEPO_cell{1,1,2,1} =reshape(block_l, [d,d,1,1,d1,1]);%right
            obj.PEPO_cell{2,1,1,1} =reshape(block_r, [d,d,d1,1,1,1]);%left
            
            obj.current_max_index=1;
            
             %initialize random boundary matrices
           
            
            
            obj.current_max_index=2;
            %obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz,d^4] ;
            %obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert,d^4] ;
            
            
            %%%%%%%%%%%%%%%%  determine 1--|--1 

            %initial guess for fsolve
            %x0 = obj.get_middle_part( {[1,2],[],[2,3],[]},[1,2,3] );
            %use contract_partial

            
           %obj.boundary_matrix_x{2,2} = rand(d1,d1);
           obj.PEPO_cell{2,1,2,1} =  obj.get_middle_part(...
                 {[1,2],[],[2,3],[]},[1,2,3],1); 

             
            if obj.testing ==1
               err = obj.calculate_error(1:3,obj.numopts) 
            end
             

            %solve 1--|--2--|--1 blocks
            part = obj.get_middle_part(...
            {[1,2],[],[3,4],[]},[1,2,3,4],0);
 
            [U,S,V] = svd( reshape(part,d1*d^2,d^2*d1));

            d2_extra = d^2;

            d2 = d^4+d2_extra;

            [U,S,V] = expand_svd(U,S,V, d2-d1*d^2 ); %still ok


            sqrt_S = diag(diag(S).^0.5);

            block_l = permute( reshape(U*sqrt_S, [d1,d,d,d2]), [2,3,1,4]);
            block_r = permute( reshape(sqrt_S*V', [d2,d,d,d1]), [2,3,1,4]); 

            obj.PEPO_cell{2,1,3,1} =reshape(block_l, [d,d,d1,1,d2,1]);%right
            obj.PEPO_cell{3,1,2,1} =reshape(block_r, [d,d,d2,1,d1,1]);%left

            if obj.testing ==1
                obj.calculate_error( [1 2 3 4],obj.numopts) 
            end   
            
            %%%%%%%

             %make 2--|--2
             obj.PEPO_cell{3,1,3,1} =  obj.get_middle_part(...
                 {[1,2,3],[],[3,4,5],[]},[1,2,3,4,5],1); 
           
            if obj.testing ==1
               err = obj.calculate_error(1:5,obj.numopts) ;
            end
             
            
            
            
            %setup new numerical prob

             %%%%%%%5solve for M2
            obj.boundary_matrix_x{1,1} = eye(1);
            obj.boundary_matrix_x{2,2} = eye(d1); 
            obj.boundary_matrix_x{3,3} = eye(d2);
            
             %%%%%%%%%%%%%5
             
             
             %prob1
            [map_1,b_map_1] = PEPO.create_map([1 2 3],obj.cycleopts);

            target_site_1 = zeros([d^2,d^2] );
            
            con_cells_1 = {...
                    {{[2,0,1,0],[1,0,2,0]},[2,1,2]},...
                    {{[2,0,2,0],[2,0,2,0]},[2,2,2]},...
                };
            
            % prob 2
            [map_2,b_map_2] = PEPO.create_map([1 2 3 4],obj.numopts);

            target2 = obj.H_exp(map_2,obj.nf);
            target_site_2 = reshape( permute(target2, site_ordering_permute(map_2.N) ), [d^2,d^2,d^2,d^2]  );       

            con_cells_2 = get_valid_contractions(obj,b_map_2, struct('max_index', obj.current_max_index));


            % prob3

            [map_3,b_map_3] = PEPO.create_map([1 2 3 4 5],obj.numopts);

            target_3 = obj.H_exp(map_3,obj.nf);
            target_site_3 = reshape( permute(target_3, site_ordering_permute(map_3.N) ), [d^2,d^2,d^2,d^2,d^2]  );       

            con_cells_3 = get_valid_contractions(obj,b_map_3, struct('max_index', obj.current_max_index));

           
            %start solver


            options = optimoptions('fsolve','Display','iter-detailed',...
                                            'Algorithm','levenberg-marquardt',...
                                            'MaxIterations',2000,...
                                            'SpecifyObjectiveGradient',true,... 
                                            'FunctionTolerance',1e-25,...
                                            'StepTolerance',1e-16,...
                                            'PlotFcn','optimplotfirstorderopt');
            x22 = obj.PEPO_cell{3,1,3,1};
            x12 = obj.PEPO_cell{2,1,3,1};
            x21 = obj.PEPO_cell{3,1,2,1};
            
            x_sizes = { size(x22), size(x12),size(x21)};   
            patterns = { [2,0,2,0], [1,0,2,0], [2,0,1,0] }; 
            begin_vec = [reshape( x22 ,1,[]) ,reshape( x12,1,[]),reshape( x21,1,[]) ];

          

            con_cells = { con_cells_1, con_cells_2 ,con_cells_3 };
            targets = {target_site_1,target_site_2,target_site_3};
            maps = { map_1,b_map_2,b_map_3 };

            [con_cells2, targets2] = optimize_con_cells(obj,maps,con_cells, patterns,targets );
           
            
           
            if obj.testing == 1
                [f0,g0] = obj.get_value_and_grad(maps,con_cells,patterns,targets,begin_vec,x_sizes);
                [f1,g1] = obj.get_value_and_grad(maps,con_cells2,patterns,targets2,begin_vec,x_sizes);
            end
            

%    
            x = fsolve( @(x) obj.get_value_and_grad(maps,con_cells2,patterns,targets2,x,x_sizes),  begin_vec , options );

            x_cell = split_x (x,x_sizes);

            x22= x_cell{1};
            x12= x_cell{2};
            x21= x_cell{3};
            
           

            obj.PEPO_cell{3,1,3,1} = x22;
            obj.PEPO_cell{2,1,3,1} = x12;
            obj.PEPO_cell{3,1,2,1} = x21;


            if obj.testing ==1
                %cyclic error improved
                err = obj.calculate_error( 1:2,obj.numopts) 
                err = obj.calculate_error( 1:3,obj.numopts) 
                err = obj.calculate_error( 1:4,obj.numopts) 

                err = obj.calculate_error(1:4,obj.cycleopts)
                err = obj.calculate_error(1:3,obj.cycleopts)
            end
             
             
             
            
 
            %%%%%%%%%%%

            obj = obj.cell2matrix() ; %save matrix form
           
        end
        


        function obj = makePEPO1d(obj)
            d = obj.dim;
            
            %todo do this in code
            obj.virtual_level_sizes_horiz = [d^0,d^2];
            obj.virtual_level_sizes_vert = [d^0,d^2];

            
            %%%%%%%%%%single site
            O_0000 =expm( 0*obj.H_1_tensor ); % eye(d);%expm( 0*obj.H_1_tensor );
            %obj.nf = trace(O_0000);
            
            obj.PEPO_cell{1,1,1,1} = reshape(  O_0000/obj.nf , [d,d,1,1,1,1] ) ;

            
            
            %%%%%%%%%%%%%% 0--|--1--|--0 and all other veriants
            obj.current_max_index = 0;
            
            
            part = obj.get_middle_part(...
                {[],[],[],[]},[1,2],0);


            [U,S,V] = svd(  reshape(part,d^2,d^2) );
            
            n_extra = 2;
            d1 = d^2+n_extra;
            
            [U,S,V] = expand_svd(U,S,V,n_extra);
            
            sqrt_S = diag(diag(S).^0.5);
            
            block_l = permute( reshape(U*sqrt_S, [1,d,d,d1]), [2,3,1,4]);
            block_r = permute( reshape(sqrt_S*V', [d1,d,d,1]), [2,3,1,4]); 
            

            obj.PEPO_cell{1,1,2,1} =reshape(block_l, [d,d,1,1,d1,1]);%right
            obj.PEPO_cell{2,1,1,1} =reshape(block_r, [d,d,d1,1,1,1]);%left
            
            obj.current_max_index=1;
            
             %initialize random boundary matrices
           
            
            
            obj.current_max_index=2;
            %obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz,d^4] ;
            %obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert,d^4] ;
            
            
            %%%%%%%%%%%%%%%%  determine 1--|--1 

            %initial guess for fsolve
            %x0 = obj.get_middle_part( {[1,2],[],[2,3],[]},[1,2,3] );
            %use contract_partial

            
           obj.boundary_matrix_x{2,2} = rand(d1,d1);
           obj.PEPO_cell{2,1,2,1} =  obj.get_middle_part(...
                 {[1,2],[],[2,3],[]},[1,2,3],1); 



            % optimize toghether

            %prob1
            [map_1,b_map_1] = PEPO.create_map([1 2 3 4],obj.cycleopts);

            target_1 = obj.H_exp(map_1,obj.nf);
            target_site_1 = reshape( permute(target_1, site_ordering_permute(map_1.N) ), [d^2,d^2,d^2]  );       

            con_cells_1 = get_valid_contractions(obj,b_map_1, struct('max_index', obj.current_max_index));

            % prob 2
            [map_2,b_map_2] = PEPO.create_map([1 2 3],obj.cycleopts);

            target2 = obj.H_exp(map_2,obj.nf);
            target_site_2 = reshape( permute(target2, site_ordering_permute(map_2.N) ), [d^2,d^2]  );       

            con_cells_2 = get_valid_contractions(obj,b_map_2, struct('max_index', obj.current_max_index));


            % prob3

            [map_3,b_map_3] = PEPO.create_map([1 2 3],obj.numopts);

            target_3 = obj.H_exp(map_3,obj.nf);
            target_site_3 = reshape( permute(target_3, site_ordering_permute(map_3.N) ), [d^2,d^2,d^2]  );       

            con_cells_3 = get_valid_contractions(obj,b_map_3, struct('max_index', obj.current_max_index));

           
            %start solver


            options = optimoptions('fsolve','Display','iter-detailed',...
                                            'Algorithm','levenberg-marquardt',...
                                            'MaxIterations',2000,...
                                            'SpecifyObjectiveGradient',true,... 
                                            'FunctionTolerance',1e-25,...
                                            'StepTolerance',1e-12);
            x11 = obj.PEPO_cell{2,1,2,1};
            M1 = obj.boundary_matrix_x{2,2};
            
            x_sizes = {size(M1)};   
            patterns = {  [1,1] }; 
            begin_vec = [reshape( M1,1,[]) ];

            con_cells = { con_cells_1 };
            targets = {target_site_1};
            maps = { b_map_1 };

% 
%             x_sizes = { size(x11), size(M1)};   
%             patterns = { [1,0,1,0], [1,1] }; 
%             begin_vec = [reshape( x11 ,1,[]) ,reshape( M1,1,[]) ];
%             
            %con_cells = { con_cells_1, con_cells_2 ,con_cells_3 };
            %targets = {target_site_1,target_site_2,target_site_3};
            %maps = { b_map_1,b_map_2,b_map_3 };

            
            [con_cells2, targets2] = optimize_con_cells(obj,maps,con_cells, patterns,targets );
            
            x = fsolve( @(x) obj.get_value_and_grad_old(maps,con_cells2,patterns,targets2,x,x_sizes),  begin_vec , options );

            x_cell = split_x (x,x_sizes);

            x1= x_cell{1};
            M1= x_cell{2};

            obj.PEPO_cell{2,1,2,1} = x1;
            obj.boundary_matrix_x{2,2} = M1;

            [U,T] = schur(M1);


            if obj.testing ==1
                %cyclic error improved
                err = obj.calculate_error( 1:2,obj.numopts) 
                err = obj.calculate_error( 1:3,obj.numopts) 
                err = obj.calculate_error( 1:4,obj.numopts) 

                err = obj.calculate_error(1:4,obj.cycleopts)
                err = obj.calculate_error(1:3,obj.cycleopts)
            end
            %%%%%%%%%end of block 1

            %solve 1--|--2--|--1 blocks
            part = obj.get_middle_part(...
            {[1,2],[],[3,4],[]},[1,2,3,4],0);
 
            [U,S,V] = svd( reshape(part,d1*d^2,d^2*d1));

            d2_extra = 4;

            d2 = d^4+d2_extra;

            [U,S,V] = expand_svd(U,S,V, d2-d1*d^2 ); %still ok


            sqrt_S = diag(diag(S).^0.5);

            block_l = permute( reshape(U*sqrt_S, [d1,d,d,d2]), [2,3,1,4]);
            block_r = permute( reshape(sqrt_S*V', [d2,d,d,d1]), [2,3,1,4]); 

            obj.PEPO_cell{2,1,3,1} =reshape(block_l, [d,d,d1,1,d2,1]);%right
            obj.PEPO_cell{3,1,2,1} =reshape(block_r, [d,d,d2,1,d1,1]);%left

            if obj.testing ==1
                obj.calculate_error( [1 2 3 4],obj.numopts) 
            end   
            
            %setup new numerical prob
            
            
            obj.PEPO_cell{3,1,3,1} =  obj.get_middle_part(...
                 {[1,2,3],[],[3,4,5],[]},[1,2,3,4,5],1); 
           
            if obj.testing ==1
               err = obj.calculate_error(1:5,obj.numopts) ;
            end
             
             %%%%%%%5solve for M2
             
            [map_2,b_map_2] = PEPO.create_map([1 2 3 4 5],obj.cycleopts);
             
            obj.boundary_matrix_x{3,3} = rand(d2,d2); 
            
            fixed_bonds_M0 = PEPO.get_fixed_bonds(b_map_2,{{[4,5],0}
                                                           {[5,1],0}});
            fixed_bonds_M1 = PEPO.get_fixed_bonds(b_map_2,{{[4,5],1}
                                                           {[5,1],1}});
            fixed_bonds_M2 = PEPO.get_fixed_bonds(b_map_2,{{[4,5],2}
                                                           {[5,1],2}});
            
            M2_target = obj.H_exp(map_2,obj.nf)...
                -obj.contract_network(b_map_2, struct('fixed',fixed_bonds_M0))...
                -obj.contract_network(b_map_2, struct('fixed',fixed_bonds_M1));
            M2_cons = obj.get_valid_contractions(b_map_2, struct('fixed',fixed_bonds_M2));
            M2_target_site = reshape( permute(M2_target, site_ordering_permute(map_2.N) ), [d^8,1]  );  
            
            [A,~]=contract_partial(obj,5, b_map_2, M2_cons );
                   
            M2 = reshape(A, [d^8, d2^2 ]) \ M2_target_site;
            
            obj.boundary_matrix_x{3,3} = reshape(M2, [d2,d2]);
            
            if obj.testing ==1
               err = obj.calculate_error(1:5,obj.cycleopts) ;
            end
            
            %%%%%%
            
            
            
             %%%%%%%%%%%%%5
             
             
             %prob1
            [map_1,b_map_1] = PEPO.create_map([1 2 3 4 5 6],obj.cycleopts);

            target_1 = obj.H_exp(map_1,obj.nf);
            target_site_1 = reshape( permute(target_1, site_ordering_permute(map_1.N) ), [d^2,d^2,d^2,d^2,d^2]  );       

            con_cells_1 = get_valid_contractions(obj,b_map_1, struct('max_index', obj.current_max_index));

            % prob 2
            [map_2,b_map_2] = PEPO.create_map([1 2 3 4 5],obj.cycleopts);

            target2 = obj.H_exp(map_2,obj.nf);
            target_site_2 = reshape( permute(target2, site_ordering_permute(map_2.N) ), [d^2,d^2,d^2,d^2]  );       

            con_cells_2 = get_valid_contractions(obj,b_map_2, struct('max_index', obj.current_max_index));


            % prob3

            [map_3,b_map_3] = PEPO.create_map([1 2 3 4 5],obj.numopts);

            target_3 = obj.H_exp(map_3,obj.nf);
            target_site_3 = reshape( permute(target_3, site_ordering_permute(map_3.N) ), [d^2,d^2,d^2,d^2,d^2]  );       

            con_cells_3 = get_valid_contractions(obj,b_map_3, struct('max_index', obj.current_max_index));

           
            %start solver


            options = optimoptions('fsolve','Display','iter-detailed',...
                                            'Algorithm','levenberg-marquardt',...
                                            'MaxIterations',2000,...
                                            'SpecifyObjectiveGradient',true,... 
                                            'FunctionTolerance',1e-25,...
                                            'StepTolerance',1e-16,...
                                            'PlotFcn','optimplotfirstorderopt');
            x22 = obj.PEPO_cell{3,1,3,1};
            M2 = obj.boundary_matrix_x{3,3};
            
            x_sizes = { size(x22), size(M2)};   
            patterns = { [2,0,2,0], [2,2] }; 
            begin_vec = [reshape( x22 ,1,[]) ,reshape( M2,1,[]) ];

            pset = 0;
            
            if pset ==0 
                con_cells = { con_cells_1, con_cells_2 ,con_cells_3 };
                targets = {target_site_1,target_site_2,target_site_3};
                maps = { b_map_1,b_map_2,b_map_3 };

                [con_cells2, targets2] = optimize_con_cells(obj,maps,con_cells, patterns,targets );
            else
                con_cells = { con_cells_2 ,con_cells_3 };
                targets = {target_site_2,target_site_3};
                maps = { b_map_2,b_map_3 };

                [con_cells2, targets2] = optimize_con_cells(obj,maps,con_cells, patterns,targets );
            end
            
           
            if obj.testing == 1
                [f0,g0] = obj.get_value_and_grad(maps,con_cells,patterns,targets,begin_vec,x_sizes);
                [f1,g1] = obj.get_value_and_grad(maps,con_cells2,patterns,targets2,begin_vec,x_sizes);
            end
            

%    
            x = fsolve( @(x) obj.get_value_and_grad(maps,con_cells2,patterns,targets2,x,x_sizes),  begin_vec , options );

            x_cell = split_x (x,x_sizes);

            x22= x_cell{1};
            M2= x_cell{2};

            obj.PEPO_cell{3,1,3,1} = x22;
            obj.boundary_matrix_x{3,3} = M2;

            [U,T] = schur(M2);


            if obj.testing ==1
                %cyclic error improved
                err = obj.calculate_error( 1:2,obj.numopts) 
                err = obj.calculate_error( 1:3,obj.numopts) 
                err = obj.calculate_error( 1:4,obj.numopts) 

                err = obj.calculate_error(1:4,obj.cycleopts)
                err = obj.calculate_error(1:3,obj.cycleopts)
            end
             
             
             
            
 
            %%%%%%%%%%%

            obj = obj.cell2matrix() ; %save matrix form
           
        end
        
        function H = H_exp(obj,map,prefactor)
            
            if nargin<3
                prefactor=1;
            end
            

            %pos_map
            % array with ones and zeros where contraction mpo are
            %should be conected
            d = obj.dim;
            H_1 = obj.H_1_tensor;
            H_2 = obj.H_2_tensor;% +... 
            %       0.5* reshape( ncon( {H_1,eye(d)}, {[-1,-3],[-2,-4]}), [d,d,d,d])+...
            %       0.5* reshape( ncon( {eye(d),H_1}, {[-1,-3],[-2,-4]}), [d,d,d,d]);

            
            H = zeros( dimension_vector( d,2*map.N));
            
            Itensor = reshape(eye(d), [d,d,1,1,1,1]);
            
            %first do all single site contributions
            for l = 1:map.N
                
                tensor_list = cell(1,map.N);
                tensor_list(:) = {Itensor};
                tensor_list{l} = H_1;
                
                leg_list_copy = map.leg_list(1,:);
                for i = 1:map.N
                    leg_list_copy{i} = leg_list_copy{i}(1:2);
                end
                
                H = H + ncon( tensor_list, leg_list_copy);
            end
            
            %do all horizontal H12
            
            for ii = 1:map.num_h_bonds
            
                arr = map.h_bonds{ii};
                n1 = arr(1);
                n2 = arr(2);
                
                if n1~=n2
                     leg_list_copy = map.leg_list(1,:);

                    index_list_n1 = map.leg_list{n1};
                    index_list_n2 = map.leg_list{n2};

                    leg_list_copy( max(n2,n1)) = []; %remove element

                    %do ij
                    new_list = [0,0,0,0];

                    new_list( [1,3] )  = index_list_n1([1,2]) ;
                    new_list( [2,4] )  = index_list_n2([1,2]) ;

                    for s=1:map.N-1
                        leg_list_copy{s} = leg_list_copy{s}(1:2);
                    end

                    leg_list_copy{ min(n1,n2) } = new_list;

                    tensor_list = cell(1,map.N-1);
                    tensor_list(:) = {Itensor};
                    tensor_list{  min(n1,n2) } = H_2;


                    H=H+ncon( tensor_list,leg_list_copy );
                end
            end
            
            %do all vertical H2
             for ii = 1:map.num_v_bonds
            
                arr = map.v_bonds{ii};
                n1 = arr(1);
                n2 = arr(2);
                
                if n1~=n2
                    leg_list_copy = map.leg_list(1,:);

                    index_list_n1 = map.leg_list{n1};
                    index_list_n2 = map.leg_list{n2};


                    leg_list_copy( max(n1,n2) ) =[]; %remove element 


                    for s=1:map.N-1
                        leg_list_copy{s} = leg_list_copy{s}(1:2);
                    end

                    %do ij
                    new_list = [0,0,0,0];


                    new_list( [1,3] )  = index_list_n1([1,2]) ;
                    new_list( [2,4] )  = index_list_n2([1,2]) ;



                    leg_list_copy{  min(n1,n2)  } = new_list;

                    tensor_list = cell(1,map.N-1);
                    tensor_list(:) = {Itensor};
                    tensor_list{   min(n1,n2) } = H_2;


                    H=H+ncon( tensor_list,leg_list_copy );
                end
                
            end
            
            H_matrix = reshape(H, [d^(map.N),d^(map.N)]);
            
            H_matrix = H_matrix -  eye(d^(map.N))*map.N*log(prefactor);
            
            H_expo = expm(H_matrix);
            
            
            H = reshape(H_expo,dimension_vector(d,2*map.N)  );
        end
        
        function contraction_cell = get_valid_contractions(obj,map, opts)
            p = inputParser;
            addParameter(p,'max_index',obj.max_index)
            addParameter(p,'matrix',0)
            addParameter(p,'fixed',zeros(map.internal_legs,1)-1)
            parse(p,opts)
            
            fixed_mask = p.Results.fixed ~=-1;
            num_fixed = sum(fixed_mask);
            

            d=obj.dim;
            
            contraction_cell = cell(0,1);
            contraction_cell_counter = 0;
            
            
            function [contraction_cell,contraction_cell_counter,vect,correct_index_set] = con_tensors(contraction_cell,contraction_cell_counter,base_level,n,stopind)

                tensor_list_indices = cell(1,map.N);

                gen_vect=  encode_index_array(n, map.internal_legs - num_fixed ,stopind)+base_level;
                vect = p.Results.fixed;
                vect(~fixed_mask) = gen_vect;
                
                correct_index_set=1;


                for i =1:map.N
                    
                    if map.is_x_border(i)||map.is_y_border(i)
                        legs = [0,0];
                        for j = 1:2
                            leg_num=map.leg_list{i}(j);
                            if leg_num > 0 
                                legs(j) = vect(leg_num);  
                            end
                        end
                        
                        if map.is_x_border(i)
                            O=isempty(obj.boundary_matrix_x{legs(1)+1,legs(2)+1});
                        else
                            O=isempty(obj.boundary_matrix_y{legs(1)+1,legs(2)+1});
                        end
                        
                        
                   
                    else
                    
                        legs=[0,0,0,0];
                        for j = 1:4
                            leg_num=map.leg_list{i}(j+2);
                            if leg_num > 0 
                                legs(j) = vect(leg_num);  
                            end
                        end
                        
                        O = isempty(obj.PEPO_cell{legs(1)+1,legs(2)+1,legs(3)+1,legs(4)+1});
                    end  
                    
                    if O==1
                        if obj.visualise==1
                            fprintf("incorrect index set \n");
                        end

                        correct_index_set=0;
                        break;
                    end

                    tensor_list_indices{i} = legs;

                end

                if correct_index_set
                    
                    contraction_cell_counter = contraction_cell_counter+1;
                    contraction_cell{contraction_cell_counter}= {tensor_list_indices,vect};
                    
                end 
            end 

            c_index = obj.cycle_index;
            if c_index == Inf

                for n=0: (p.Results.max_index+1)^(map.internal_legs-num_fixed)  -1

                      [contraction_cell,contraction_cell_counter,~] = con_tensors(contraction_cell,contraction_cell_counter,0,n,p.Results.max_index);
                end

            else
                fprintf("new");
                for n=0: (obj.cycle_index )^(map.internal_legs-num_fixed)-1
                      [contraction_cell,contraction_cell_counter,vect, correct_index_set] = con_tensors(contraction_cell,contraction_cell_counter,0,n,obj.cycle_index-1);                      
                end

                fprintf("cycle_ind");
                for n=0: ( (p.Results.max_index-obj.cycle_index) +1)^(map.internal_legs-num_fixed)-1
                      [contraction_cell,contraction_cell_counter,vect,correct_index_set] = con_tensors(contraction_cell,contraction_cell_counter,c_index,n,p.Results.max_index-obj.cycle_index);      
                 end

            end

        end
         
        function tensors = fetch_PEPO_cells(obj,map,legs,patterns,xs)
            
            if nargin <=3
               patterns = [];
            end
            
            num_patterns = size(patterns,2);
         
            
            tensors=cell(1,map.N);
            for n=1:map.N
                leg = legs{n};

                matched_pattern = 0;
                
                for ii=1:num_patterns
                    if same_pattern(leg,patterns{ii})==1
                        tensors{n} = xs{ii};
                        matched_pattern = 1;
                        break;
                    end
                end
                
                if matched_pattern == 0
                    if map.is_x_border(n)
                       tensors{n} = obj.boundary_matrix_x{leg(1)+1,leg(2)+1};
                    elseif map.is_y_border(n)
                       tensors{n} = obj.boundary_matrix_y{leg(1)+1,leg(2)+1};
                    else
                        tensors{n} = obj.PEPO_cell{leg(1)+1,leg(2)+1,leg(3)+1,leg(4)+1};
                    end
                end

            end
        end

        function [A,x0_shape]= contract_partial(obj,num, map, con_cells , x, pattern )  

            %xsizes: how to split x into different components 
            %patterns: match x with pattern given pattern for subs

            
            % create new contraction list without num and legs connected to
            % num as last indices
            l=map.h_bond_l_lookup{num};
            r=map.h_bond_r_lookup{num};
            u=map.v_bond_u_lookup{num};
            d=map.v_bond_d_lookup{num};


            con_list_cpy = map.leg_list;

            ii = 0;

            if ~isempty(l)
                ii = ii+1; 
                pair = map.h_bonds{l};
                other = pair(1);
                
                if map.is_x_border(other)
                    con_list_cpy{other}(2) = -(map.external_legs+ii);
                else
                    con_list_cpy{other}(2+3) = -(map.external_legs+ii);
                end

            end

            if ~isempty(u)
                ii = ii+1; 
                pair = map.v_bonds{u};
                other = pair(1);
                if map.is_y_border(other)
                    con_list_cpy{other}(2) = -(map.external_legs+ii);
                else
                    con_list_cpy{other}(2+4) = -(map.external_legs+ii);
                end
            end

            if ~isempty(r)
                ii = ii+1; 
                pair = map.h_bonds{r};
                other = pair(2);
                if map.is_x_border(other)
                    con_list_cpy{other}(1) = -(map.external_legs+ii);
                else
                    con_list_cpy{other}(2+1) = -(map.external_legs+ii);
                end
            end

            if ~isempty(d)
                ii = ii+1; 
                pair = map.v_bonds{d};
                other = pair(2);
                if map.is_y_border(other)
                    con_list_cpy{other}(1) = -(map.external_legs+ii);
                else
                    con_list_cpy{other}(2+2) = -(map.external_legs+ii);
                end
            end

            num_legs_cpy = con_list_cpy{num};

            x0_list = con_list_cpy{ num};
            con_list_cpy( num) = [];


            final_order = -1:-1: -(map.external_legs+ii);
            seq = 1:map.internal_legs; %contraction sequence


            final_external_missing =sort( - num_legs_cpy(num_legs_cpy<0),'descend');
            final_internal_missing = sort( num_legs_cpy(num_legs_cpy>0),'descend');

            for l = 1:size(final_external_missing,2)
                final_order(  final_external_missing(l)  )= [];
            end

            for l = 1:size(final_internal_missing,2)
                seq(  final_internal_missing(l)  )= [];
            end
            
            
            

            % do actual contractions
            A = [];
            
            
            for con_cell_index =1:size(con_cells,2)
                legs = con_cells{con_cell_index}{1};
                
                if nargin <5
                    temp_list = fetch_PEPO_cells(obj,map,legs);
                else
                    temp_list = fetch_PEPO_cells(obj,map,legs,pattern,x);
                end
                
                x0 = temp_list{num}; 
                x0_shape = size(x0);
                
                temp_list(num) = []; %remove x0
                
                if isempty(A)
                    A = ncon( temp_list,con_list_cpy,seq,final_order  );
                else
                    
                    A=A+ncon( temp_list,con_list_cpy,seq,final_order  );
                end
                
            end

            
            
            %put in site ordering
            if map.is_x_border(num) || map.is_x_border(num)
                perm_vector = [site_ordering_permute(map.N2); ((2*map.N2+1):size(size(A),2)).']; 
                A = reshape( permute(A, perm_vector)  , [],prod(x0_shape(1:2)));
            else
                a_size = size(A);

                ext_legs = x0_list(3:end);
                smallest_ind = min(-ext_legs( ext_legs<0 ));
                
                idx = find(final_order == -( smallest_ind -1));
                
                num1=2*(map.N2-1);
                num2=idx;
                num3 = size(size(A),2)-ii;
                
                size1 = a_size(1:num1);%ij indices 
                size2 = a_size(num1+1:num2); %external legs before
                size3 = a_size(num2+1:num3); %external legs after
                size4 = a_size(num3+1:end); %bond to x0
                
                perm_vector = [site_ordering_permute(map.N2-1); ((2*map.N2-1):size(size(A),2)).']; 
                A = reshape( permute(A, perm_vector)  ,[prod(size1),prod(size2),1,prod(size3),prod(size4) ]);
            end

        end
        
        function [F,G] = get_value_and_grad(obj,maps,con_cells_cell,patterns,targets,x0,x_sizes)
            
            num_sub_probs = size(maps,2);
            
            F_cell= cell(num_sub_probs,1);
            G_cell = cell(num_sub_probs,1);
            
            
            F_total_size = 0;
            G_total_size = numel(x0);
            
            nargout_val = nargout;
            
            for sub_prob = 1:num_sub_probs
                con_cells = con_cells_cell{sub_prob};
                target = targets{sub_prob};
                map = maps{sub_prob};
                
                F_total_size = F_total_size+numel(target);
                
                
                num_x = size(x_sizes,2);

                %x= reshape(x,[],1);

                x_cell = cell(num_x,1);


                curr=0;
                for i =1:num_x
                    num_elem = prod(x_sizes{i});
                    x_cell{i} = reshape(x0(curr+1:curr+num_elem),x_sizes{i});
                    curr = curr+num_elem;
                end

                total_g_params = curr;
                
                F_sub_buffer = cell(1,size(con_cells,2)  );
                G_sub_buffer = cell(1,size(con_cells,2)  );
                

                %could be parforred
                parfor con_cell_index =1:size(con_cells,2)


                    legs = con_cells{con_cell_index}{1};
                    temp_list_1 = fetch_PEPO_cells(obj,map,legs,patterns,x_cell);

                    A1  = ncon( temp_list_1,map.leg_list);
                    %F_sub = F_sub + reshape(  permute(A1,site_ordering_permute(map.N2)), size(target));

                    
                    perm_vect = [site_ordering_permute(map.N2), 2*map.N2+1:  size(size(A1),2) ] ;
                    
                    F_sub_buffer{con_cell_index} = reshape( permute(A1,perm_vect), size(target));
                    
                    if nargout_val  > 1 %calculate gradient

                        
                        
                        G_sub_cell = cell(num_x,1);
                        for i =1:num_x
                            num_elem = prod(x_sizes{i});
                            G_sub_cell{i} = zeros(numel(target),num_elem  );
                        end
                        
                        
                        for pat_num = 1:num_x

                            x=x_cell{pat_num};

                            if size(patterns{pat_num},2)==2 %boundary matrix

                                grad_total_size = [ numel(target), numel(x) ];
                                Grad_total = zeros( grad_total_size );

                                for ii = 1:size(legs,2)
                                    if  same_pattern(legs{ii}, patterns{ pat_num} )
                                        [Ai,~]=contract_partial(obj,ii, map, con_cells(con_cell_index) , x_cell,patterns  ) ;
                                        Grad_total = Grad_total + Ai;
                                    end
                                end

                            else
                                size_x_red = size(x,[1,2,3,4,5,6]);
                                d2 = obj.dim^2;
                                size_x_red(1) = d2;
                                size_x_red(2)=[];


                                grad_total_size = [size(target),size_x_red];
                                Grad_total = zeros( grad_total_size );

                                
                                

                                num = 0;

                                for ii = 1:size(legs,2)
                                    if  same_pattern(legs{ii}, patterns{ pat_num} )

                                        index_before = num;
                                        index_after = map.N2-num-1;

                                        
                                        external_sizes = numel(target)/d2^(map.N2);
                                        
                                        
                                        [Ai,~]=contract_partial(obj,ii, map, con_cells(con_cell_index) , x_cell,patterns  ) ;

                                        size_x = size_x_red(2:5);
                                        non_connected = map.leg_list{ii}<0;
                                        non_connected = non_connected(3:end);
                                        size_x_internal = size_x;
                                        size_x_internal( non_connected) = 1;
                                        
                                        size_x_external = size_x;
                                        size_x_external(~non_connected) = 1;
                                        
                                        
                                        
                                        ai_size= size(Ai);
                                        external_size = ai_size(2:4);
                                        size_curr = external_sizes/prod(external_size);
                                        
                                        Grad_total = reshape(Grad_total, [d2^index_before, d2,d2^index_after,external_size(1),size_curr,external_size(3), size_x_red]     );

                                        
                                        
                                        
                                        Ai_res = reshape(Ai, [d2^index_before,1,d2^index_after,external_size,1, size_x_internal  ] );
                                                      
                                          
                                        
                                       
                                        if isequal([0,0,0,0] ,non_connected)
                                            for iii =1:size_x_red(1)
                                                Grad_total(:,iii,:, :,1,:, iii,:,:,:,:) = Grad_total(:,iii,:, :,1,:, iii,:,:,:,:)...
                                                    + Ai_res(:,1,:, :,1,:,  1,:,:,:,:);
                                            end
                                        elseif isequal([0,1,1,1] ,non_connected)
                                            for iii =1:size_x_red(1)
                                                for iiii =1:size_curr
                                                    [i1,i2,i3,i4] = ind2sub(size_x_external,iiii);
                                                    
                                                    Grad_total(  :,iii,:, :,iiii,:, iii,:,i2,i3,i4) = Grad_total(:,iii,:, :,iiii,:, iii,:,i2,i3,i4)...
                                                        + Ai_res(:,1,  :, :,1,   :, 1,  :,1,1,1);
                                                end
                                            end
                                        elseif  isequal([1,1,0,1] ,non_connected)
                                            for iii =1:size_x_red(1)
                                                for iiii =1:size_curr
                                                    [i1,i2,i3,i4] = ind2sub(size_x_external,iiii);
                                                    Grad_total(:,iii,:, :,iiii,:, iii,i1,i2,:,i4) = Grad_total(:,iii,:, :,iiii,:, iii,i1,i2,:,i4)...
                                                        + Ai_res(:,1,:, :,1,:,   1,   1,1,:,1);
                                                end
                                            end
                                        elseif isequal([0,1,0,1] ,non_connected)
                                            for iii =1:size_x_red(1)
                                                for iiii =1:size_curr
                                                    [i1,i2,i3,i4] = ind2sub(size_x_external,iiii);
                                                    
                                                    Grad_total(  :,iii,:, :,iiii,:, iii,:,i2,:,i4) = Grad_total(:,iii,:, :,iiii,:, iii,:,i2,:,i4)...
                                                        + Ai_res(:,1,  :, :,1,   :, 1,  :,1, :,1);
                                                end
                                            end    
                                        else
                                            non_connected
                                            error("not implemented")

                                         end


                                        Grad_total =  reshape(Grad_total, grad_total_size );
                                    
                                    end

                                    if size(legs{ii},2)==4 %not boundary matrix
                                        num = num+1;
                                    end
                                end
                            end

                            G_sub_cell{pat_num} = G_sub_cell{pat_num} +reshape( Grad_total, size(G_sub_cell{pat_num}));
                        end

                        G_sub_local = zeros( [ numel(target) ,total_g_params ]);

                        curr=0;
                        for i =1:num_x
                            num_elem = prod(x_sizes{i});
                            G_sub_local(:,curr+1:curr+num_elem) = G_sub_cell{i}  ;
                            curr = curr+num_elem;

                        end
                        
                        G_sub_buffer{con_cell_index} = G_sub_local;

                    end
                end
                
                
                
                F_sub = -target;
                if nargout_val > 1 
                    G_sub = 0;
                end
                
                
                for con_cell_index =1:size(con_cells,2)
                      F_sub = F_sub+  F_sub_buffer{con_cell_index};
                      if nargout_val > 1 
                            G_sub = G_sub+  G_sub_buffer{con_cell_index};
                      end
                end
                
                %put back toghether
                F_cell{sub_prob} = F_sub;
                if nargout_val > 1 
                    G_cell{sub_prob} = G_sub;
                end
            
            end
            
           F = zeros(F_total_size,1);
           if nargout_val  > 1 
              G = zeros(F_total_size,G_total_size);
           end


           curr = 0;

           for sub_prob = 1:num_sub_probs
                num_elem =  numel(F_cell{sub_prob})  ;
                F(curr+1:curr+num_elem) = reshape(F_cell{sub_prob},[num_elem,1]) ;

                if nargout > 1 
                    G(curr+1:curr+num_elem,:) = G_cell{sub_prob};
                end

                curr = curr+num_elem;

           end

            
        end

        function [con_cells_cell2, targets] = optimize_con_cells(obj,maps,con_cells_cell, patterns,targets )
            
            num_sub_probs = size(maps,2);
            

            con_cells_cell2 = cell(1,num_sub_probs);
            
            
            for sub_prob = 1:num_sub_probs
                con_cells = con_cells_cell{sub_prob};
                con_cells2 = cell(1,1);
                con_cells2_counter = 1;
                
                target = targets{sub_prob};
                map = maps{sub_prob};
               
                num_x = size(patterns,2);

                has_matched_pattern =0;
                
                
                for con_cell_index =1:size(con_cells,2)
                    legs = con_cells{con_cell_index}{1};
                    for pat_num = 1:num_x
                        for ii = 1:size(legs,2)
                            if  same_pattern(legs{ii}, patterns{ pat_num} )
                                has_matched_pattern=1;
                                break;
                            end
                        end
                        if has_matched_pattern
                            break;
                        end
                    end
                    
                    if has_matched_pattern %keep in new list
                        con_cells2{con_cells2_counter} = con_cells{con_cell_index};
                        con_cells2_counter =  con_cells2_counter+ 1;
                    else %remove and change target
                        temp_list_1 = fetch_PEPO_cells(obj,map,legs);
                        
                        A1  = ncon( temp_list_1,map.leg_list);
                        target = target - reshape(  permute(A1,site_ordering_permute(map.N2)), size(target));
                        
                    end
                end
                
                con_cells_cell2{1,sub_prob} = con_cells2;
                targets{sub_prob} = target;
                
            end
        end
        
        
        function M = contract_network(obj,map, opts)
            %generate all index sets for a given configuration
                        
%             struct('max_index', [] ,'matrix', [] ,'fixed', []);
            

            p = inputParser;
            addParameter(p,'max_index',obj.max_index)
            addParameter(p,'matrix',0)
            addParameter(p,'fixed',zeros(map.internal_legs,1)-1)
            parse(p,opts)
            
            
            M = zeros( dimension_vector(obj.dim,2*map.N2) );
            
            if p.Results.matrix == 0 
                correct_index_sets = get_valid_contractions(obj,map, opts);
                
                for i = 1:size(correct_index_sets,2)
                    iset = correct_index_sets{i};
                    vect = iset{2};
                    legs = iset{1};
                    
                    tensors = obj.fetch_PEPO_cells(map,legs);

                    M = M + ncon( tensors, map.leg_list );
                    
                end
                

            else
                tensor_list = cell(1,map.N);
                for i = 1:map.N
                    T = obj.PEPO_matrix;
                    connections = map.leg_list{i};
                    %only keep sublevel 0 for the given tensors
                    if connections(1+2)<0
                       T = T(:,:,1,:,:,:); 
                    end
                    if connections(2+2)<0
                       T = T(:,:,:,1,:,:); 
                    end
                    if connections(3+2)<0
                       T = T(:,:,:,:,1,:); 
                    end
                    if connections(4+2)<0
                       T = T(:,:,:,:,:,1); 
                    end
                    tensor_list{i} = T;
                end
                
                M = ncon_optim( tensor_list, map.leg_list  );
                
            end
            
        end
        
        function obj = cell2matrix(obj)
                        
            d = obj.dim;
            
            size_arr_horiz = obj.virtual_level_sizes_horiz;

            start_index_H = zeros(obj.max_index+2, 1);
            start_index_H(1) = 1;
            ind = 1;
            for i = 2:obj.max_index + 2
                ind = ind + size_arr_horiz(i-1);
                start_index_H(i) = ind;
            end


            totaldimensionH = start_index_H(end) - 1;

            
            function y = getH(i)
               y=  start_index_H(i):start_index_H(i + 1)-1;
            end
            
            size_arr_vert = obj.virtual_level_sizes_vert;

            start_index_V = zeros(obj.max_index+2, 1);
            start_index_V(1) = 1;
            ind = 1;
            for i = 2:obj.max_index + 2
                ind = ind + size_arr_vert(i-1);
                start_index_V(i) = ind;
            end


            totaldimensionV = start_index_V(end) - 1;

            
            function y = getV(i)
               y=  start_index_V(i):start_index_V(i + 1)-1;
            end
            
            

            T = zeros(d,d,totaldimensionH,totaldimensionV,totaldimensionH,totaldimensionV);
            
            %sparsem = ndSparse(sparse(d^2, totaldimensionH^2*totaldimensionV^2));
            %T = reshape(sparsem, [d,d,totaldimensionH,totaldimensionV,totaldimensionH,totaldimensionV]);
            
            
            %move all existing tensors to matrix
            for i1 = 1:obj.max_index + 1
                for i2 = 1:obj.max_index + 1
                    for i3 = 1:obj.max_index + 1
                        for i4 = 1:obj.max_index + 1
                            %trace the spins for the environment
                            cell =  obj.PEPO_cell{i1,i2,i3,i4};
                            if length(cell)~=0 
                                T(:,:,getH(i1), getV(i2), getH(i3),getV(i4)) = obj.PEPO_cell{i1,i2,i3,i4}; %ncon(  { obj.PEPO_cell{i1,i2,i3,i4} } , {[1,1,-1,-2,-3,-4]} );
                            end
                        end
                    end
                end
            end

           
            obj.PEPO_matrix = T;
            

            %K=reshape(T(:,1,1,:),[45,45]);


%             obj.left = zeros(1, totaldimension);
%             obj.left(1) = 1;
%             obj.right = zeros(totaldimension, 1);
%             obj.right(1) = 1;



            
        end
        
        function get_reflection(obj)
            a=ncon( {obj.PEPO_cell{1,1,1,2}},{[1,1,-1,-2,-3,-4]})-ncon( {obj.PEPO_cell{1,2,1,1}},{[1,1,-1,-4,-3,-2]});
            b=ncon( {obj.PEPO_cell{1,1,2,1}},{[1,1,-1,-2,-3,-4]})-ncon( {obj.PEPO_cell{2,1,1,1}},{[1,1,-3,-2,-1,-4]});

            disp(a);
            
        end
        
        function [err,prefact] = calculate_error(obj,nummap,opts)
            [map,b_map ] = PEPO.create_map(nummap,opts);
            
            d=obj.dim;
            

            H_matrix=H_exp(obj,map,obj.nf);
            
            Contraction = obj.contract_network(b_map,struct('max_index', obj.max_index))  ;
          
            
            b = reshape(  H_matrix, [ d^(map.N2), d^(map.N2)]);
            a = reshape( Contraction, [ d^(map.N2), d^(map.N2)]);
           
            p = 2;

            [~, S1, ~] = svds(a-b, 30);


            sum_1 = (sum(diag(S1).^p))^(1 / p);


            [~, S2, ~] = svds(b, 30);


            sum_2 = (sum(diag(S2).^p))^(1 / p);
            
           
            prefact = obj.nf^map.N*( sum_2);
            
            % for real values, multiply both with obj.nf^(map.N)*trace_a =
            % prefact^N

            err =  sum_1/ sum_2   ;
        end
        
        function [A,B,G1,lambda1] = vumps(obj,chimax)
            

            %todo check these params
            opts.charges='regular';
            opts.dynamical='off';
            opts.dyncharges=0;
            opts.schmidtcut=1e-10;
            opts.chimax=350;
            %opts.disp='iter';
            opts.disp='none';
            opts.tolmax=1e-4; %1e-4
            opts.tolfactor=1e4;
            opts.minit=1;
            opts.dyniter=5;
            opts.truncate=0;
            opts.method='vumps';
            opts.save=0;

            %opts.method = 'qr';
         

            opts.plot='on';
            opts.maxit=1000;
            opts.tolfixed=1e-12;
            
            %put into vumps format
            

            T = obj.PEPO_matrix;
         
            hdim = size(T,3);
            vdim = size(T,4);

            %upper vumps zipper
            M = ncon(  {T },  {[1,1,-1,-2,-3,-4] } );
            %lower vumps zipper
            

            
            o.legs=4;
            o.group='none';
            o.dims = size( M ) ;
            o.var = M;
            
          
            

            O.type = 'mpo';
            O.mpo=o;
            

            [A,G1,lambda1,~,~]= Vumps(O,chimax,[],opts);
            
            %correct estimate for inversion sym?
            
            GL = G1{1};GR = G1{2};Ac = A{4}; 
            
            m.legs=4;
            m.group='none';
            m.dims = size( M ) ;
            m.var = M;
            
            
            opts.krylovdim=100; opts.tol=1e-14;
            %opts.disp='iter-detailed'; %opts.reorth='force';

            function x = vumps_under(x)
                x=TensorContract({x,GL,m,GR},{[1,2,5],[1,3,-1],[-2,4,2,3],[-3,4,5]});
            end
            
            [B,lambda2,err2]=TensorEigs(@(x) vumps_under(x), TensorConj(Ac),1,'lm',opts);
      
            
            
            
            %[B_l,C_l,~]=TensorDecRight(Bc,'polar');
            
            %[B_r,C_r,~]=TensorDecLeft(Bc,'polar');
            

            
%             accopts.method = 'qr';
%             accopts.tol = 1e-16;
%             
%             [B,err]=VumpsSolveACC(B2,B_c,accopts);
            
            
            
        end
        
        function part = get_middle_part(obj,inv_maps,map,perm)
             %calculates residual hamiltonian and inverts the legs given in
             %inv_maps in the order {left,up,right,down}. left and up
             %should be increasing indices to central pepo cell, right and
             %down decreasing.
            
             if nargin<4
                perm = 1;
             end
             
            d=obj.dim;
           
            
            
            [map,b_map] = PEPO.create_map(map, obj.numopts);
            
            Tensor = obj.H_exp(map,obj.nf)-obj.contract_network(map,struct('max_index', obj.current_max_index) );

            Tensor_site = permute(Tensor, site_ordering_permute(map.N));
            
            vector_sizes = [0,0,0,0,0];
            for l=1:2
                i_map = inv_maps{l};
                if size(i_map,1) ~= 0
                    vector_sizes(l) = sum( i_map~=0  )-1; 
                end 
            end
            for l=3:4
                i_map = inv_maps{l};
                if size(i_map,1)~= 0
                    vector_sizes(l+1) = sum( i_map~=0  )-1; 
                end 
            end
            vector_sizes(3) = map.N-sum(vector_sizes);
            vector_sizes = d.^(2* vector_sizes );
            
            
            
            Tensor_site = reshape(Tensor_site, vector_sizes  );
            
            if obj.testing
                Tensor_site_cpy = Tensor_site;
            end

            
            function [l_chain_site,bond_size] = get_l_chain(l_map,l_map_2) 
                l_map_2 = PEPO.create_map(l_map_2,obj.numopts);
                l_map = PEPO.create_map(l_map, obj.numopts);
                
                l_tensors = cell(1,l_map.N-1);
            
                last_index = 0;
                for i = 1:l_map.N-1
                    index_set = [0,0,0,0];

                    index_set_mask = l_map.leg_list{i}(3:end)>0; 
                    non_empty = find(index_set_mask);

                    if size(non_empty,2)==1
                        last_index = non_empty(1);
                        index_set(last_index) = i;
                    else
                        prev_index = mod(last_index+2,4);
                        if non_empty(1)==prev_index
                            index_set(non_empty(1))=i-1;
                            index_set(non_empty(2))=i;
                            last_index = non_empty(2);
                        else
                            index_set(non_empty(2))=i-1;
                            index_set(non_empty(1))=i;
                            last_index = non_empty(1);
                        end
                    end
                
                    l_tensors{i} = obj.PEPO_cell{index_set(1)+1,index_set(2)+1,index_set(3)+1,index_set(4)+1};
                end

                leg_list_copy = l_map_2.leg_list;
                l_chain =  ncon( l_tensors,  leg_list_copy   );

                l_chain_size = size(l_chain);
                l_size = d^(2*l_map_2.N);
                
                bond_size= numel(l_chain)/ l_size ;
                l_chain_size_new = [l_chain_size(1:2*l_map_2.N),bond_size] ;
                
                l_chain = reshape(l_chain, l_chain_size_new );
                
                
                l_permute = [ site_ordering_permute(l_map_2.N); 2*l_map_2.N+1];        
                l_chain_site = reshape( permute( l_chain, l_permute ), [l_size,bond_size] );
            end
            
            function [r_chain_site,bond_size] = get_r_chain(r_map,r_map_2) 

                r_map_2 = PEPO.create_map(r_map_2,obj.numopts);
                r_map = PEPO.create_map(r_map,obj.numopts);
                
                r_tensors = cell(1,r_map.N-1);
            
                last_index = 0;
                for i = r_map.N:-1:2
                    index_set = [0,0,0,0];

                    index_set_mask = r_map.leg_list{i}(3:end)>0; 
                    non_empty = find(index_set_mask);

                    if size(non_empty,2)==1
                        last_index = non_empty(1);
                        index_set(last_index) = r_map.N- i+1;
                    else
                        prev_index = mod(last_index+2,4);
                        if non_empty(1)==prev_index
                            index_set(non_empty(1))=r_map.N-(i);
                            index_set(non_empty(2))=r_map.N-(i-1);
                            last_index = non_empty(2);
                        else
                            index_set(non_empty(2))=r_map.N-(i);
                            index_set(non_empty(1))=r_map.N-(i-1);
                            last_index = non_empty(1);
                        end
                    end

                    r_tensors{i-1} = obj.PEPO_cell{index_set(1)+1,index_set(2)+1,index_set(3)+1,index_set(4)+1};
                end

               
                leg_list_copy = r_map_2.leg_list;
                r_chain = ncon( r_tensors,  leg_list_copy   );

                r_chain_size = size(r_chain);
                r_size = d^(2*r_map_2.N);
                
                bond_size= numel(r_chain)/ r_size ;
                
                r_chain_new_size = [r_chain_size(1:2*r_map_2.N) ,bond_size]; 

                r_chain = reshape(r_chain, r_chain_new_size  );
                
                r_permute = [ 2*r_map_2.N+1;  site_ordering_permute(r_map_2.N)] ;               
                r_chain_site = reshape( permute( r_chain, r_permute ), [bond_size,r_size] ); 
            end
            
            function [x1,x2] = renumber(x,central,renumber)
               mask =  x~=0 ;
               
               mask_central = x==central;
               
               x_size = size(x,1);
               y = reshape( x(mask), [],1 );

               [~,y] = sort(y);
               y = reshape(y,x_size,[]);
               
               x1=x;
               x1(mask) = y;
               x2 = x1;
               x2(mask_central) = 0;
               
               if renumber==1
                    x2(~mask_central) = x2(~mask_central)-1;
               end

            end
            
            if obj.testing ==1
               tensorarr = cell(1,4);
               for ii = 1:4
                   tensorarr{ii} = [1];
               end
            end
            
            for leg=1:4
                inv_map = inv_maps{leg};
                if size(inv_map,1) ~= 0
                    
                    if leg<3
                        central = max(inv_map);
                        
                        [i_map,i_map_2] = renumber(inv_map,central,0);
                        [l_chain_site,bond_size] = get_l_chain(i_map,i_map_2);
                        
                        if obj.testing ==1
                           tensorarr{leg} = l_chain_site; 
                        end
                        
                        if leg ==1 
                            Tensor_site = reshape(Tensor_site, vector_sizes(1),[]);
                            Tensor_site2 = l_chain_site\Tensor_site;
                            vector_sizes(1)=bond_size;
                            Tensor_site = reshape(Tensor_site2,vector_sizes);
                        else
                            Tensor_site = reshape( permute(Tensor_site, [2,1,3,4,5] ) , vector_sizes(2),[]);
                            Tensor_site = l_chain_site\Tensor_site;
                            vector_sizes(2)=bond_size;
                            Tensor_site = permute( reshape(Tensor_site,[vector_sizes(2),vector_sizes(1),vector_sizes(3),vector_sizes(4),vector_sizes(5)]), [2,1,3,4,5]);
                        end 
                    else
                        central = min(inv_map(inv_map~=0));
                        

                        
                        [i_map,i_map_2] = renumber(inv_map,central,1);
                        [r_chain_site,bond_size] = get_r_chain(i_map,i_map_2);                        
                        if obj.testing ==1
                           tensorarr{leg} = r_chain_site; 
                        end
                        if leg ==4
                            Tensor_site = reshape(Tensor_site, [],vector_sizes(5));
                            Tensor_site = Tensor_site/r_chain_site;
                            vector_sizes(5)=bond_size;
                            Tensor_site = reshape(Tensor_site,vector_sizes);
                        else
                            Tensor_site = reshape( permute(Tensor_site, [1,2,3,5,4] ) , [],vector_sizes(4));
                            Tensor_site2 = Tensor_site/r_chain_site;
                            vector_sizes(4)=bond_size;
                            Tensor_site = permute( reshape(Tensor_site2, [vector_sizes(1),vector_sizes(2),vector_sizes(3),vector_sizes(5),vector_sizes(4)]), [1,2,3,5,4]);
                        end 
                    end
            
                end
                
            end
            
            if obj.testing
                Z=ncon( {tensorarr{1},tensorarr{2},Tensor_site,tensorarr{3},tensorarr{4} },{ [-1,1],[-2,2],[1,2,-3,3,4],[3,-4],[4,-5]  }  );
            
                ZZ= Z-Tensor_site_cpy;
            end
            
            
            if perm
                part = reshape( permute( Tensor_site, [3,1,2,4,5]), [d,d,vector_sizes(1),vector_sizes(2),vector_sizes(4),vector_sizes(5)]);
            else
                part = Tensor_site;
            end
            
        end
        
        function [mag, corr_length,delta] = get_expectation (obj,X,chimax)
            [A,B,G,lambda]  = vumps(obj,chimax);
            
            
            T = obj.PEPO_matrix;
            
            M = ncon(  {T},  {[1,1,-1,-2,-3,-4]} );
            
            m.legs=4;
            m.group='none';
            m.dims = size( M ) ;
            m.var = M;
            
            O =  ncon ( {T,X}, {[1,2,-1,-2,-3,-4],[1,2]}  );
            o.legs=4;
            o.group='none';
            o.dims = size( O ) ;
            o.var = O;
            
            %transfereigs edited !
            
            GL = G{1};GR = G{2};Ac = A{4}; 
          
            
            %should be lambda?
            [x,~]=TensorContract({B,GL,Ac,m,GR},...
                    {[1,2,6],[1,3,4],[4,5,8],[5,7,2,3],[8,7,6]});
                
            [y,~]=TensorContract({B,GL,Ac,o,GR},...
                    {[1,2,6],[1,3,4],[4,5,8],[5,7,2,3],[8,7,6]});    
            
            mag = y/x;
            
            
            
            function x= transfer_up(x)
                [x,~]=TensorContract({GL,x,m,GR},...
                    {[-1,3,4],[4,5,8],[5,7,-2,3],[8,7,-3]});
            end
            
            opts.krylovdim=100; opts.tol=1e-14;
            
            [rho,f] = TensorEigs(@(x) transfer_up(x),A{4},5,'lm',opts);

            eps_i = -log(abs(f));
            corr_length = eps_i(1+1)-eps_i(1);
            
            delta = eps_i(4+1)-eps_i(2+1); %same as https://arxiv.org/pdf/1907.08603.pdf
        end
        
    end
    
    methods (Access= private, Static=true)
        
        function x = to_vumps_order(x)
           x = ncon({x},{ [-2,-3,-4,-1]}  );
        end
        
    end
    
    methods (Static=true)
        function [map,boundary_map] = create_map(pos_map,opts,internal)
            %h_cyclic -> right column is boundary matrix mx
            %V_cyclic -> under row is boundary matrix my
            %matrices on the boundary should recieve the higher number than
            %the others
            
            if nargin < 2
               opts.numbered = 0; 
               opts.v_cyclic = 0;
               opts.h_cyclic = 0;
            end
            
            if nargin < 3
               internal = 0; %recursive call
            end

            
            if ~isfield(opts,'numbered')
                opts.numbered=0;
            end
            
            if ~isfield(opts,'v_cyclic')
                opts.v_cyclic=0;
            end
          
            if ~isfield(opts,'h_cyclic')
                opts.h_cyclic=0;
            end
                
            %number the location of operators from up to down and left to
            %right, and create toghether with it a leg_list for ncon

            [m,n] = size(pos_map);
            map.m=m;
            map.n=n;
            
            map.pos_lookup = {};

            
            if opts.numbered == 0
                counter = 1;
                %boundary matrices have highest numbers
                for x =1:n-1
                    for y =1:m-1
                        if pos_map(y,x)==1
                            pos_map(y,x) = counter;
                            map.pos_lookup{counter} = [y,x];
                            counter = counter +1;
                        end
                    end
                end
                
               
                for y =1:m-1
                    if pos_map(y,n)==1
                        pos_map(y,n) = counter;
                        map.pos_lookup{counter} = [y,n];
                        counter = counter +1;
                    end
                end
                
                for x =1:n
                    if pos_map(m,x)==1
                        pos_map(m,x) = counter;
                        map.pos_lookup{counter} = [m,x];
                        counter = counter +1;
                    end
                end
                
                
                N = counter-1;
            else %todo more checking
                N = 0;
                for x =1:n
                    for y =1:m
                        if pos_map(y,x)~=0
                            counter = pos_map(y,x);
                            %pos_map(y,x) = counter;
                            map.pos_lookup{counter} = [y,x];
                            if counter > N
                               N = counter; 
                            end
                        end
                    end
                end
            end

            %N number of physical sites

            map.N=N;

            map.num_map=pos_map;

            map.h_bonds = {}; 
            map.v_bonds = {};
            
            %left right up down bonds
            map.h_bond_l_lookup = cell(N,1);
            map.h_bond_r_lookup = cell(N,1);
            map.v_bond_u_lookup = cell(N,1);
            map.v_bond_d_lookup = cell(N,1);
            
            map.is_x_border = zeros(map.N,1);
            map.is_y_border = zeros(map.N,1);

            leg_list = cell(1,N);
            leg_list(:) = { [0,0,0,0,0,0]  };

            %do the horizontal internal bonds

            internal_counter = 1;

            for num = 1:N
               coor = map.pos_lookup{num};
               x= coor(2); y = coor(1);
               
               if x==n
                   if opts.h_cyclic==1
                       next_x = 1;
                       map.is_x_border(n)=1;    
                   else
                       continue; %skip this round
                   end
               else
                  next_x = x+1; 
               end
                
               if pos_map(y,x) ~=0 && pos_map(y,next_x) ~=0
                    n1 = pos_map(y,x);
                    n2 = pos_map(y,next_x);

                    leg_list{n1}(5) = internal_counter;
                    leg_list{n2}(3) = internal_counter;

                    
                    map.h_bonds{internal_counter} = [n1,n2];
                    
                    %save bonds per site
                    map.h_bond_r_lookup{n1} = [ map.h_bond_r_lookup{n1} , internal_counter ];
                    map.h_bond_l_lookup{n2} = [ map.h_bond_l_lookup{n2} , internal_counter ];
                    
                    internal_counter=internal_counter+1;

                    
                    %fprintf( "hor %d-%d\n",pos_map(y,x), pos_map(y,x+1));
               end
            end
            
            map.num_h_bonds = internal_counter-1;
            
            %vertical internal bonds

            for num = 1:N
               coor = map.pos_lookup{num};
               x= coor(2); y = coor(1);
               
               if y==m
                   if opts.v_cyclic==1
                       next_y = 1;
                       map.is_y_border(n)=1;
                   else
                       continue; %skip this round
                   end
               else
                  next_y = y+1; 
               end
               

                if pos_map(y,x) ~=0 && pos_map(next_y,x) ~=0
                    n1 = pos_map(y,x);
                    n2 = pos_map(next_y,x);

                    leg_list{n1}(6) = internal_counter;
                    leg_list{n2}(4) = internal_counter;

                    v_number = internal_counter-map.num_h_bonds;
                    map.v_bonds{ v_number} = [n1,n2];
                    
                    %save bonds per site
                    map.v_bond_d_lookup{n1} = [ map.v_bond_d_lookup{n1} , v_number ];
                    map.v_bond_u_lookup{n2} = [ map.v_bond_u_lookup{n2} , v_number ];
                    
                    
                    internal_counter=internal_counter+1;

                    %fprintf( "vert %d-%d\n",pos_map(y,x), pos_map(y+1,x));
                end
                
            end

            
            map.num_v_bonds = internal_counter-map.num_h_bonds-1;
            map.internal_legs=internal_counter-1;

            
            
            if internal ==0
            
                total_borders = sum(map.is_x_border)+sum(map.is_y_border);
                map.N2 = map.N-total_borders;

                if sum(map.is_x_border(1:map.N2))~=0
                    error("boundary matrices should have largest numbers")
                end

                if sum(map.is_y_border(1:map.N2))~=0
                    error("boundary matrices should have largest numbers")
                end

               

            else %boundary row and column already deleted.
                map.N2 =  map.N;
                map.is_x_border= map.is_x_border*0;
                map.is_y_border= map.is_y_border*0;
            end
            
        
            %number ij according to number
            external_counter=1;
            
            for N1 =1:map.N2
                
                    leg_list{N1}(1)= -(N1)  ;
                    leg_list{N1}(2)= -(N1+ map.N2 );

                    external_counter = external_counter+1;
            end
            
            for N1 =map.N2+1:map.N     
                   leg_list{N1} =  leg_list{N1}( leg_list{N1} ~= 0);
            end

            external_counter= 2*(map.N2)+1;

            %number all other indices
            
            for i=1:N
                if map.is_x_border(i)==0 && map.is_y_border(i)==0
                    for j=3:6
                        if leg_list{i}(j)==0
                            leg_list{i}(j) = -external_counter;
                            external_counter = external_counter+1;
                        end
                    end
                end
            end

            map.external_legs=external_counter-1;
            map.leg_list= leg_list;

            
            map.map = "true";
            
            if internal == 0
                if opts.h_cyclic || opts.v_cyclic
                    map2 = PEPO.create_map(  map.num_map(1:end-opts.v_cyclic,1:end-opts.h_cyclic) ,opts,1);
                
                    boundary_map = map;
                    boundary_map.boundary_map=1;
                    
                    map = map2;
                    map.boundary_map=0;
                    
                    return;
                else
                    boundary_map = map;
                    return;%normal non cyclic map
                end
            else
                
                return;
            end

        end

        function fixed_bonds = get_fixed_bonds(map,bonds)
            %bonds is nx1 cell with 2x1 cells with {[n1,n2], virt_level} 
            
            num_fixed = size(bonds,2);
            
            
            fixed_bonds= zeros(map.internal_legs,1)-1;
            
            for i=1:num_fixed
                bond = bonds{i};
                n = bond{1};
                n1= n(1); n2=n(2);
                virt_level = bond{2};
                
                [h_bond,~]=intersect([map.h_bond_r_lookup{n1},map.h_bond_l_lookup{n1}],[map.h_bond_r_lookup{n2},map.h_bond_l_lookup{n2}]);
                [v_bond,~]=intersect([map.v_bond_u_lookup{n1},map.v_bond_d_lookup{n1}],[map.v_bond_u_lookup{n2},map.v_bond_d_lookup{n2}]);
                
                for ii = 1:size(h_bond,1)
                    fixed_bonds( h_bond(ii) ) = virt_level;
                end
                
                for ii = 1:size(v_bond,1)
                    fixed_bonds( v_bond(ii)+map.num_h_bonds ) = virt_level;
                end
                
                pot_bond_v{1} = map.v_bond_u_lookup{n1};
                pot_bond_v{2} = map.v_bond_u_lookup{n2};
                
                
            end
            
        end
    end
end


%helper functions

function p = dimension_vector(d,n,leftright)
%helper function to create a 1xn vector  [ left,d,d,..,d,right]
%if left/right are not supplied/0, this is omitted

if nargin<3
    p = zeros(1,n);
    p = p+d;
    return
    
else
    
    p = zeros(1,n+2);
    p = p+d;
    p(1) = leftright(1);
    p(end) = leftright(2);
end

end

function y = encode_index_array(n,len,max_num)
% this takes a single number and splits it into the composite indices. The
i =1;
y = zeros(len,1);

while n ~= 0
    [n,r] = readOne(n,max_num+1);
    y(i) = r;
    i = i+1;
end

end

function [s,r] = readOne(s,d)
r = rem(s,d);
s = (s-r)/d;
end

function p = site_ordering_permute(n)
% changes from |left i1 i2 ... j1 j2.. right> to |left i1 j1 i2 j2 ...
% right>
p = zeros(2*n,1);
%p(1)=1;
%p(2*n+2)=2*n+2;
for i = 1:n
    p(2*i-1)=i;
end
for i = 1:n
    p(2*i) = n+i ;
end
end

function bool = same_pattern(leg1,leg2)
    bool = 0;
    if size(leg1)==size(leg2)
        if leg1==leg2
           bool=1; 
        end
    end
end

function x_cell = split_x (x,x_sizes)

    num_x = size(x_sizes,2);
    x_cell = cell(num_x,1);

    curr=0;
    for i =1:num_x
        num_elem = prod(x_sizes{i});
        x_cell{i} = reshape(x(curr+1:curr+num_elem),x_sizes{i});
        curr = curr+num_elem;
    end
            
end


function [U2,S2,V2] = expand_svd(U,S,V,n)
    dim = size(S,1);
    U2=zeros(dim,dim+n);
    V2= zeros(dim,dim+n);
    S2=zeros(dim+n);
    
    if dim < dim+n
        U2(1:dim,1:dim) = U;
        V2(1:dim,1:dim) = V;
        S2(1:dim,1:dim) = S;
    else
        U2 = U(:,1:dim+n);
        V2 = V(:,1:dim+n);
        S2 = S(1:dim+n,1:dim+n);
        
    end
end

function C= rotate_rhs(B,num,physical_dims)
    %B = reshape(B,physical_dims);

    perm = 1:size(size(B),2);
    perm(num)=[];
    perm = [perm,num];
  
    C = permute(B, perm );
end
