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
        
        function obj = makePEPO_3(obj)
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
            sqrt_S = diag(diag(S).^0.5);
            
            block_l = permute( reshape(U*sqrt_S, [1,d,d,d^2]), [2,3,1,4]);
            block_r = permute( reshape(sqrt_S*V', [d^2,d,d,1]), [2,3,1,4]); 
            

            obj.PEPO_cell{1,1,2,1} =reshape(block_l, [d,d,1,1,d^2,1]);%right
            obj.PEPO_cell{2,1,1,1} =reshape(block_r, [d,d,d^2,1,1,1]);%left
            
            %%%%%%%%%%%%%%%%%create 0--|--1--|--1--|--0 and variants

            obj.current_max_index=1;
            
            block_11 = obj.get_middle_part( {[1,2],[],[2,3],[]},[1,2,3] );

            obj.PEPO_cell{2,1,2,1}= reshape(block_11,[d,d,d^2,1,d^2,1]); 

            if obj.testing ==1
                obj.calculate_error( PEPO.create_map([1 2 3],obj.numopts)) 
            end

                                
            %%%%%%%%%%%%% create  0--|--1--|--2--|--1--|--0

            %horizontal
            obj.current_max_index=2;

            obj.virtual_level_sizes_horiz = [obj.virtual_level_sizes_horiz,d^4] 
            obj.virtual_level_sizes_vert = [obj.virtual_level_sizes_vert,d^4] 
            

            part = obj.get_middle_part(...
                {[1,2],[],[3,4],[]},[1,2,3,4],0);

            [U,S,V] = svd( reshape(part,d^2*d^2,d^2*d^2) );
            sqrt_S = diag(diag(S).^0.5);

            block_l = permute( reshape(U*sqrt_S, [d^2,d,d,d^4]), [2,3,1,4]);
            block_r = permute( reshape(sqrt_S*V', [d^4,d,d,d^2]), [2,3,1,4]); 

            obj.PEPO_cell{2,1,3,1}= reshape(block_l,[d,d,d^2,1,d^4,1]);

            obj.PEPO_cell{3,1,2,1}= reshape(block_r,[d,d,d^4,1,d^2,1]);


            if obj.testing ==1
                obj.calculate_error( PEPO.create_map([1 2 3 4],obj.numopts)) 
            end
                    
            converged =0;
            
            while converged ~=1
                [map,b_map] = PEPO.create_map([1 2 3 4],obj.cycleopts);

                fixed_bonds = PEPO.get_fixed_bonds(b_map,{ {[3,4],1},...
                                                           {[4,1],1} });

                target = obj.H_exp(map,obj.nf)-...
                    obj.contract_network(b_map,struct('max_index', obj.current_max_index, 'fixed', fixed_bonds  ))
                
                
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
        
         
        function M = contract_network(obj,map, opts)
            %generate all index sets for a given configuration
                        
%             struct('max_index', [] ,'matrix', [] ,'fixed', []);
            
            p = inputParser;
            addParameter(p,'max_index',obj.max_index)
            addParameter(p,'matrix',0)
            addParameter(p,'fixed',zeros(map.internal_legs,1)-1)
            parse(p,opts)
            
          
            fixed_mask = p.Results.fixed ~=-1;
            num_fixed = sum(fixed_mask);
            
            if obj.visualise ==1
                new_map = zeros( 2*map.m+1, 2*map.n+1 )+NaN;
                for j=1:map.N
                   coor= map.pos_lookup{j} ;
                   new_y = 2+ 2*(coor(1)-1);
                   new_x = 2+ 2*(coor(2)-1);
                   
                   new_map( new_y, new_x)= -j;
                end
            end
            
            d=obj.dim;
            M = zeros( dimension_vector( d,2*map.N) );
            
            function [M,vect,correct_index_set] = con_tensors(M,base_level,n,stopind)

                tensor_list = cell(1,map.N);

                gen_vect=  encode_index_array(n, map.internal_legs - num_fixed ,stopind)+base_level;
                vect = p.Results.fixed;
                vect(~fixed_mask) = gen_vect;
                
                

                correct_index_set=1;

                if obj.visualise==1
                    disp(vect);
                    map_copy = new_map(:,:);
                end

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
                            O=obj.boundary_matrix_x{legs(1)+1,legs(2)+1};
                        else
                            O=obj.boundary_matrix_y{legs(1)+1,legs(2)+1};
                        end
                        
                        
                   
                    else
                    
                        legs=[0,0,0,0];
                        for j = 1:4
                            leg_num=map.leg_list{i}(j+2);
                            if leg_num > 0 
                                legs(j) = vect(leg_num);  
                            end
                        end
                        
                        O = obj.PEPO_cell{legs(1)+1,legs(2)+1,legs(3)+1,legs(4)+1};
                    end  


                    
                    
                    if length(O)==0
                        if obj.visualise==1
                            fprintf("incorrect index set \n");
                        end

                        correct_index_set=0;
                        break;
                    end

                    if ~correct_index_set

                        break;
                    end

                    if obj.visualise==1 %probably doens't work for cyclic
                       coor= map.pos_lookup{i} ;
                       new_y = 2+ 2*(coor(1)-1);
                       new_x = 2+ 2*(coor(2)-1);

                       map_copy(new_y,new_x-1) = legs(1);
                       map_copy(new_y,new_x+1) = legs(3);
                       map_copy(new_y-1,new_x) = legs(2);
                       map_copy(new_y+1,new_x) = legs(4);

                    end
                    
                    tensor_list{i} = O;

                end

                if correct_index_set
                    if obj.visualise==1
                        disp(map_copy);
                    end

                    M=M+ncon_optim( tensor_list, map.leg_list);
                end 
            end 

            if p.Results.matrix == 0 
                c_index = obj.cycle_index;
                if c_index == Inf
                  
                    for n=0: (p.Results.max_index+1)^(map.internal_legs-num_fixed)  -1
                
                          [M,~] = con_tensors(M,0,n,p.Results.max_index);
                    end

                else
                    fprintf("new");
                    for n=0: (obj.cycle_index )^(map.internal_legs-num_fixed)-1
                          [M,vect, correct_index_set] = con_tensors(M,0,n,obj.cycle_index-1);                      
                    end

                    fprintf("cycle_ind");
                    for n=0: ( (p.Results.max_index-obj.cycle_index) +1)^(map.internal_legs-num_fixed)-1
                          [M,vect,correct_index_set] = con_tensors(M,c_index,n,p.Results.max_index-obj.cycle_index);      
                     end
                    
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
        
        function [err,prefact] = calculate_error(obj,map,matrix)
            if nargin <3 
               matrix = 0; 
            end
            
            d=obj.dim;
            
            

            H_matrix=H_exp(obj,map,obj.nf);
            
            Contraction = obj.contract_network(map,struct('max_index', obj.max_index,'matrix',matrix))  ;
            %Contraction=obj.contract_network(map_arg,obj.max_index,1);

            
            
            b = reshape(  H_matrix, [ d^(map.N), d^(map.N)]);
            a = reshape( Contraction, [ d^(map.N), d^(map.N)]);
            
            %eig_a= eigs(a,4);
            
            
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
           
            
            
            map = PEPO.create_map(map, obj.numopts);
            
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

            
            function [l_chain_site,l_size] = get_l_chain(l_map,l_map_2) 
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
                l_chain_new_size = [l_chain_size(1:2*l_map_2.N) ,l_size]; 

                l_chain = reshape(l_chain, l_chain_new_size  );
                
                
                l_permute = [ site_ordering_permute(l_map_2.N); 2*l_map_2.N+1];        
                l_chain_site = reshape( permute( l_chain, l_permute ), [l_size,l_size] );
            end
            
            function [r_chain_site,r_size] = get_r_chain(r_map,r_map_2) 

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
                r_chain_new_size = [r_chain_size(1:2*r_map_2.N) ,r_size]; 

                r_chain = reshape(r_chain, r_chain_new_size  );
                
                r_permute = [ 2*r_map_2.N+1;  site_ordering_permute(r_map_2.N)] ;               
                r_chain_site = reshape( permute( r_chain, r_permute ), [r_size,r_size] ); 
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
                        [l_chain_site,l_size] = get_l_chain(i_map,i_map_2);
                        
                        if obj.testing ==1
                           tensorarr{leg} = l_chain_site; 
                        end
                        
                        if leg ==1 
                            Tensor_site = reshape(Tensor_site, vector_sizes(1),[]);
                            Tensor_site2 = l_chain_site\Tensor_site;
                            Tensor_site = reshape(Tensor_site2,vector_sizes);
                        else
                            Tensor_site = reshape( permute(Tensor_site, [2,1,3,4,5] ) , vector_sizes(2),[]);
                            Tensor_site = l_chain_site\Tensor_site;
                            Tensor_site = permute( reshape(Tensor_site,[vector_sizes(2),vector_sizes(1),vector_sizes(3),vector_sizes(4),vector_sizes(5)]), [2,1,3,4,5]);
                        end 
                    else
                        central = min(inv_map(inv_map~=0));
                        

                        
                        [i_map,i_map_2] = renumber(inv_map,central,1);
                        [r_chain_site,r_size] = get_r_chain(i_map,i_map_2);                        
                        if obj.testing ==1
                           tensorarr{leg} = r_chain_site; 
                        end
                        if leg ==4
                            Tensor_site = reshape(Tensor_site, [],vector_sizes(5));
                            Tensor_site = Tensor_site/r_chain_site;
                            Tensor_site = reshape(Tensor_site,vector_sizes);
                        else
                            Tensor_site = reshape( permute(Tensor_site, [1,2,3,5,4] ) , [],vector_sizes(4));
                            Tensor_site2 = Tensor_site/r_chain_site;
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
                    leg_list{N1}(2)= -(N1+ map.N );

                    external_counter = external_counter+1;
            end
            
            for N1 =map.N2+1:map.N     
                   leg_list{N1} =  leg_list{N1}( leg_list{N1} ~= 0);
            end

            external_counter= 2*(map.N)+1;

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