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
    end
    
    
        
    methods
        function obj = PEPO(d,H_1_tensor,H_2_tensor,max_index,type,opts)
            obj.dim = d;
            obj.H_1_tensor = H_1_tensor;
            obj.H_2_tensor = H_2_tensor;
            obj.PEPO_cell = cell( max_index+1,max_index+1,max_index+1,max_index+1 );
            obj.type = type;
            obj.max_index = max_index;
            
            p = inputParser;
            addParameter(p,'testing',0)
            addParameter(p,'visualise',0)
            parse(p,opts)
            
            obj.testing= p.Results.testing;
            obj.visualise = p.Results.visualise;
            
            
            obj = obj.makePEPO();
        end
        
        function obj = makePEPO(obj)
            d = obj.dim;
            
            %todo do this in code
            obj.virtual_level_sizes_horiz = [d^0,d^2,d^2];
            obj.virtual_level_sizes_vert = [d^0,d^2,d^4];
            %obj.virtual_level_sizes_horiz = [d^0,d^2];
            %obj.virtual_level_sizes_vert = [d^0,d^2];
            
            %%%%%%%%%%single site
            O_0000 = expm( obj.H_1_tensor );
            obj.nf = trace(O_0000);
            
            obj.PEPO_cell{1,1,1,1} = reshape(  O_0000/obj.nf , [d,d,1,1,1,1] ) ;

            %%%%%%%%%%%%%% 0--|--1--|--0 and all other veriants
            obj.current_max_index = 0;
            
            map = PEPO.create_map([1 1]);

            Tensor_010 = obj.H_exp(map,obj.nf) -...
                obj.contract_network(map,obj.current_max_index);

            Tensor_010_site = reshape(  permute(Tensor_010, site_ordering_permute(map.N)),...
                            [d^2,d^2] );

            [U,S,V] = svd( Tensor_010_site);

            sqrt_S = diag(diag(S).^0.5);
            
            block_01 = reshape(U*sqrt_S, [d,d,d^2]);
            block_10 = permute( reshape(sqrt_S*V', [d^2,d,d]), [2,3,1]); 

            obj.PEPO_cell{1,1,2,1} =reshape(block_01, [d,d,1,1,d^2,1]);%right
            obj.PEPO_cell{1,1,1,2} =reshape(block_01, [d,d,1,1,1,d^2]);%down

            obj.PEPO_cell{2,1,1,1} =reshape(block_10, [d,d,d^2,1,1,1]);%left
            obj.PEPO_cell{1,2,1,1} =reshape(block_10, [d,d,1,d^2,1,1]);%up


            if obj.testing==1
                err = ncon( { obj.PEPO_cell{1,1,2,1},  obj.PEPO_cell{2,1,1,1}  },{ [-1,-3,-5,-6,1,-7],[-2,-4,1,-8,-9,-10]  }  )-Tensor_010;
                fprintf(  "err 010 %.4e \n",  eigs(  reshape(err, [d^map.N,d^map.N]),1) ); 
            end
            %end
            
            %%%%%%%%%%%%%%%%%create 0--|--1--|--1--|--0 and variants
            obj.current_max_index=1;
            
            part = obj.get_middle_part( [1,2],[2,3],[1,2,3] );
            
            block_11 = permute( reshape(part, [d^2,d,d,d^2]), [2,3,1,4]);

            obj.PEPO_cell{2,1,2,1}= reshape(block_11,[d,d,d^2,1,d^2,1]);
            obj.PEPO_cell{2,1,1,2}= reshape(block_11,[d,d,d^2,1,1,d^2]);
            obj.PEPO_cell{1,2,2,1}= reshape(block_11,[d,d,1,d^2,d^2,1]);
            obj.PEPO_cell{1,2,1,2}= reshape(block_11,[d,d,1,d^2,1,d^2]);

            %this is between 2 01 blocks

            part = obj.get_middle_part( [1,2],[2;
                                               3],[0 3; 
                                                   1 2]);

            block_11 = permute( reshape(part, [d^2,d,d,d^2]), [2,3,1,4]);
            obj.PEPO_cell{2,2,1,1}= reshape(block_11,[d,d,d^2,d^2,1,1]);

            %this is between 2 10 blocks

            part = obj.get_middle_part( [1;
                                         2],[2,3],[2 3; 
                                                   1 0]); 
            block_11 = permute( reshape(part, [d^2,d,d,d^2]), [2,3,1,4]);

            obj.PEPO_cell{1,1,2,2}= reshape(block_11,[d,d,1,1,d^2,d^2]);

            
            
            obj.max_index=1;



%             %  correct for one loop
             
            map = PEPO.create_map([1 1;
                                   1 1]);
            
            Tensor_1111 = obj.H_exp(map,obj.nf)-...
                obj.contract_network(map,obj.current_max_index);
            
%             err = eigs(reshape(Tensor_1111,[d^(map.N),d^(map.N)]),1);   
%             fprintf("block err %.4e",abs(err));
%             
%             if abs(err) > 1e-8
                Tensor_1111_site = reshape(  permute(Tensor_1111, site_ordering_permute(map.N)),...
                                [d^4,d^4] );

                [U,S,V] = svd(Tensor_1111_site );
                sqrt_S = S.^0.5;

                l = reshape( permute(reshape(U*sqrt_S, [d^2,d^2,d^2,d^2]), [3,1,2,4]), [d^4,d^4]);
                [U2,S2,V2] = svd(l);

                sqrt_S2 = S2.^0.5;

                lu = reshape(  U2*sqrt_S2, [d^2,d,d,d^4]);
                ld =  reshape(  sqrt_S2*V2', [d^4,d,d,d^2]);

                r = reshape( permute(reshape(sqrt_S*V', [d^2,d^2,d^2,d^2]), [1,3,4,2]), [d^4,d^4]);
                [U3,S3,V3] = svd(r);

                sqrt_S3 = S3.^0.5;

                ru = reshape( U3*sqrt_S3, [d^2,d,d,d^4]);
                rd =  reshape(  sqrt_S3*V3', [d^4,d,d,d^2]);

                lu_2 = permute(lu,[2,3,1,4]);
                ld_2 = permute(ld,[2,3,1,4]);
                ru_2 = permute(ru,[2,3,1,4]);
                rd_2 = permute(rd,[2,3,4,1]);

                obj.PEPO_cell{1,1,3,3}= reshape(lu_2,[d,d,1,1,d^2,d^4]);
                obj.PEPO_cell{1,3,3,1}= reshape(ld_2,[d,d,1,d^4,d^2,1]);
                obj.PEPO_cell{3,1,1,3}= reshape(ru_2,[d,d,d^2,1,1,d^4]);
                obj.PEPO_cell{3,3,1,1}= reshape(rd_2,[d,d,d^2,d^4,1,1]);

                if obj.testing
                    %err1 = eigs(reshape(ncon(  {lu,ld,ru,rd}, { [1,-1,-5,2],[2,-2,-6,3],[1,-3,-7,4],[4,-4,-8,3] })-Tensor_1111,[d^4,d^4]),1);
                    err = eigs(reshape(....
                        ncon( {obj.PEPO_cell{1,1,3,3},obj.PEPO_cell{1,3,3,1},obj.PEPO_cell{3,1,1,3},obj.PEPO_cell{3,3,1,1}},...               
                             map.leg_list )-Tensor_1111...
                          ,[d^4,d^4]),1);
                      fprintf("err decomposing block %.4e",err);
                end
%                      
                
                obj = obj.cell2matrix() ; %save matrix form
            %end            
        end
        
        function H = H_exp(obj,map,prefactor)
            
            if nargin<3
                prefactor=1;
            end
            
            %map either contains pos_map or (map from parse_map function)   
            map = PEPO.parse_map(map);

            %pos_map
            % array with ones and zeros where contraction mpo are
            %should be conected
            d = obj.dim;
            H_1 = obj.H_1_tensor;
            H_2 = obj.H_2_tensor;
            
            
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
  
            
            for x =1:map.n-1
                for y =1:map.m
                    if map.num_map(y,x) ~=0 && map.num_map(y,x+1) ~=0
                        a = map.num_map(y,x);
                        b = map.num_map(y,x+1);
                        
                        n1 = min(a,b);
                        n2 = max(a,b);
                        
                        
                        leg_list_copy = map.leg_list(1,:);
                        
                        index_list_n1 = map.leg_list{n1};
                        index_list_n2 = map.leg_list{n2};
                        
                        leg_list_copy(n2) = []; %remove element
                        
                        %do ij
                        new_list = [0,0,0,0];
                        
                        new_list( [1,3] )  = index_list_n1([1,2]) ;
                        new_list( [2,4] )  = index_list_n2([1,2]) ;
                        
                        for s=1:map.N-1
                            leg_list_copy{s} = leg_list_copy{s}(1:2);
                        end
                        
                        leg_list_copy{n1} = new_list;
                        
                        tensor_list = cell(1,map.N-1);
                        tensor_list(:) = {Itensor};
                        tensor_list{n1} = H_2;
                        
                        
                        H=H+ncon( tensor_list,leg_list_copy );
                    end
                end
            end
            
            %do all vertical H2
            for x =1:map.n
                for y =1:map.m-1
                    if map.num_map(y,x) ~=0 && map.num_map(y+1,x) ~=0
                        
                        a = map.num_map(y,x);
                        b = map.num_map(y+1,x);
                        
                        n1 = min(a,b);
                        n2 = max(a,b);
                        

                        leg_list_copy = map.leg_list(1,:);
                        
                        index_list_n1 = map.leg_list{n1};
                        index_list_n2 = map.leg_list{n2};
                        
                        
                        leg_list_copy( n2) =[]; %remove element 
 
                        
                        for s=1:map.N-1
                            leg_list_copy{s} = leg_list_copy{s}(1:2);
                        end
                        
                        %do ij
                        new_list = [0,0,0,0];
                        
                        
                        new_list( [1,3] )  = index_list_n1([1,2]) ;
                        new_list( [2,4] )  = index_list_n2([1,2]) ;

                            
                        
                        leg_list_copy{  n1  } = new_list;
                        
                        tensor_list = cell(1,map.N-1);
                        tensor_list(:) = {Itensor};
                        tensor_list{n1} = H_2;
                        
                        
                        H=H+ncon( tensor_list,leg_list_copy );
                    end
                end
            end
            
            H_matrix = reshape(H, [d^(map.N),d^(map.N)]);
            
            H_matrix = H_matrix -  eye(d^(map.N))*map.N*log(prefactor);
            
            H_expo = expm(H_matrix);
            
            
            H = reshape(H_expo,dimension_vector(d,2*map.N)  );
        end
        
        function M = contract_network(obj,map, max_index,matrix,borders)
            %generate all index sets for a given configuration
                        
            if nargin < 4 
                matrix = 0;
            end
            
            if nargin < 5
                if matrix==1
                    error("not implemented")
                end
                
                borders = [0,0];
            end
            
            
            map = PEPO.parse_map(map);
            
            if obj.visualise ==1
                new_map = zeros( 2*map.m+1, 2*map.n+1 )+NaN;
                for j=1:map.N
                   coor= map.lookup{j} ;
                   new_y = 2+ 2*(coor(1)-1);
                   new_x = 2+ 2*(coor(2)-1);
                   
                   new_map( new_y, new_x)= -j;

                end
            end
            
            d=obj.dim;
            M = zeros( dimension_vector( d,2*map.N) );
            
            
            if matrix == 0 
                %with zero included
                for n=0: (max_index+1)^map.internal_legs-1

                    %todo
                    tensor_list = cell(1,map.N);

                    vect=  encode_index_array(n, map.internal_legs,max_index);

                    correct_index_set=1;

                    if obj.visualise==1
                        disp(vect);
                        map_copy = new_map(:,:);
                    end

                    for i =1:map.N
                        legs=[0,0,0,0];
                        for j = 1:4
                            leg_num=map.leg_list{i}(j+2);
                            if leg_num > 0 
                                legs(j) = vect(leg_num);  
                            end
                        end  

                        
                        

                        O = obj.PEPO_cell{legs(1)+1,legs(2)+1,legs(3)+1,legs(4)+1};

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

                        if obj.visualise==1
                           coor= map.lookup{i} ;
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
        
        function [err,prefact] = calculate_error(obj,map)
            d=obj.dim;
            
            map = PEPO.parse_map(map);
            

            H_matrix=H_exp(obj,map,obj.nf);
            
            Contraction = obj.contract_network(map,obj.max_index,0);
            %Contraction=obj.contract_network(map_arg,obj.max_index,1);

            
            
            a = reshape(  H_matrix, [ d^(map.N), d^(map.N)]);
            b = reshape( Contraction, [ d^(map.N), d^(map.N)]);
            
            eig_a= eigs(a,1);
           
            prefact = obj.nf^map.N*(eig_a);
            
            % for real values, multiply both with obj.nf^(map.N)*trace_a =
            % prefact^N

            err =  ( eigs(  (a-b)  ,1) )/ eig_a   ;
        end
        
        function [A,G,lambda,ctr,error] = vumps(obj)
            

            %todo check these params
            opts.charges='regular';
            opts.dynamical='off';
            opts.dyncharges=0;
            opts.schmidtcut=1e-6;
            opts.chimax=350;
            opts.disp='iter';
            opts.tolmax=1e-4; %1e-4
            opts.tolfactor=1e5;
            opts.minit=1;
            opts.dyniter=5;
            opts.truncate=0;
            opts.method='vumps';
            opts.save=0;

            %opts.method = 'qr';
            

            opts.plot='on';
            opts.maxit=1000;
            opts.tolfixed=1e-8;
            
            %put into vumps format
            

            T = obj.PEPO_matrix;
         
            hdim = size(T,3);
            vdim = size(T,4);
            
            %put auxilary indices at the end for vumps and reshape to peps
            %format
            %M = ncon(  {T},  {[1,1,-1,-2,-3,-4]} );
            
            %M = reshape(ncon(  {T},  {[-1,-4,-2,-3,-5,-6]} ),...
            %                        [ obj.dim, hdim*vdim, obj.dim, hdim*vdim ]);
            
            %M = permute( reshape( T, [obj.dim,obj.dim,hdim*vdim,hdim*vdim]),...
            %    [1,3,2,4]);
            
            M = permute (T, [3,4,5,6,1,2]);
            
            %M = ncon(  {T,conj(T)},  {[1,2,-1,-3,-5,-7],[2,1,-2,-4,-6,-8]} );
            %M = reshape(M, [hdim^2,vdim^2,hdim^2,vdim^2]);
           % M2 = full(M);
            
            o.legs=4;
            o.group='none';
            o.dims = size( M ) ;
            o.var = M;

            O.type = 'mpo';
            O.mpo=o;


            [A,G,lambda,ctr,error]=Vumps(O,5,[],opts);

        end
        
        function part = get_middle_part(obj,l_map,r_map,map)
             %inverts l_map and r_map 
             %l_map ; largest index is open
             %only for chain like structures
            
            d=obj.dim;
           
            r_map_2 = r_map(2:end);
            r_map_2 = r_map_2-min(r_map_2)+1;
            r_map_2 = PEPO.create_map(r_map_2,1);
            r_map = r_map-min(r_map)+1;
            r_map = PEPO.create_map(r_map,1);
            
            l_map_2 = PEPO.create_map(l_map(1:end-1),1);
            l_map = PEPO.create_map(l_map,1);
            
            
            
            map = PEPO.create_map(map,1);
            
            % get the correct tensors for left inverse. Last one is not
            % added
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
            
            % saame for right side. First one doesn't count
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
            
            
            %%%%%%%%%%%%%todo apply inverses
            
            Tensor = obj.H_exp(map,obj.nf)-obj.contract_network(map,obj.current_max_index);

            Tensor_site = reshape(  permute(Tensor, site_ordering_permute(map.N)),...
                            l_size,[] );
            
                        
                        
            l_permute = [ site_ordering_permute(l_map_2.N); 2*l_map_2.N+1];        
            r_permute = [ 2*r_map_2.N+1;  site_ordering_permute(r_map_2.N)] ;               
             
            l_chain_site = reshape( permute( l_chain, l_permute ), [l_size,l_size] ); 
            r_chain_site = reshape( permute( r_chain, r_permute ), [r_size,r_size] ); 
         
            x = l_chain_site\Tensor_site;
            x = reshape(x, [],r_size);
            y = x/r_chain_site;
            
            part = reshape(y, l_size,[],r_size);
        end
        
        function m = get_expectation (obj,X)
            [A,G,lambda,ctr,error] = vumps(obj);
            
            M = permute( reshape( T, [obj.dim,obj.dim,hdim*vdim,hdim*vdim]),...
                [1,3,2,4]);
            
            
            
            T =   ncon ( {M,X }, {[1,-1,2,-2],[1,2]}  );
            O =  ncon ( {M}, {[1,-1,2],[1,2]}  );
            
            
            Ac = A{3};
            Gl = G{1};
            Gr = G{2};
            
            overlap1=TensorContract({G{1},A{4},G{2},TensorConj(A{4})},{[1,4:A{1}.legs+1,2],[2,3],[3,4:A{1}.legs+2],[1,A{1}.legs+2]});
            
      
            
        end
        
    end
    
    methods (Access= private, Static=true)
        function map = parse_map(map)
            if isfield(map,"pos_map")
                map = PEPO.create_map(map.pos_map);
            else
                if isfield(map,"map")
                    map=map;
                else
                    error("wrong input struct");
                end
            end
        end
        
        function x = to_vumps_order(x)
           x = ncon({x},{ [-2,-3,-4,-1]}  );
        end
        
    end
    
    methods (Static=true)
        function map = create_map(pos_map,numbered)
            if nargin < 2
               numbered = 0; 
            end
                
            
            %number the location of operators from up to down and left to
            %right, and create toghether with it a leg_list for ncon

            [m,n] = size(pos_map);
            map.m=m;
            map.n=n;

            map.lookup = {};

            
            if numbered == 0
            
                counter = 1;
                for x =1:n
                    for y =1:m
                        if pos_map(y,x)==1
                            pos_map(y,x) = counter;
                            map.lookup{counter} = [y,x];
                            counter = counter +1;
                        end
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
                            map.lookup{counter} = [y,x];
                            if counter > N
                               N = counter; 
                            end
                        end
                    end
                end
            end

            

            map.N=N;

            map.num_map=pos_map;


            leg_list = cell(1,N);
            leg_list(:) = { [0,0,0,0,0,0]  };

            %do the horizontal internal bonds

            internal_counter = 1;

            for x =1:n-1
                for y =1:m
                    if pos_map(y,x) ~=0 && pos_map(y,x+1) ~=0
                        n1 = pos_map(y,x);
                        n2 = pos_map(y,x+1);

                        leg_list{n1}(5) = internal_counter;
                        leg_list{n2}(3) = internal_counter;

                        internal_counter=internal_counter+1;

                        %fprintf( "hor %d-%d\n",pos_map(y,x), pos_map(y,x+1));
                    end
                end
            end

            %vertical internal bonds

            for x =1:n
                for y =1:m-1
                    if pos_map(y,x) ~=0 && pos_map(y+1,x) ~=0
                        n1 = pos_map(y,x);
                        n2 = pos_map(y+1,x);

                        leg_list{n1}(6) = internal_counter;
                        leg_list{n2}(4) = internal_counter;

                        internal_counter=internal_counter+1;

                        %fprintf( "vert %d-%d\n",pos_map(y,x), pos_map(y+1,x));
                    end
                end
            end

            map.internal_legs=internal_counter-1;

            %number ij according to number
            external_counter=1;

            for x =1:n
                for y =1:m
                    N1 = pos_map(y,x);
                    if N1 ~=0
                        leg_list{N1}(1)= -N1;
                        leg_list{N1}(2)= -(N1+N);

                        external_counter = external_counter+1;
                    end
                end
            end

            external_counter= 2*N+1;

            %number all other indices
            
            for i=1:N
                for j=3:6
                    if leg_list{i}(j)==0
                        leg_list{i}(j) = -external_counter;
                        external_counter = external_counter+1;
                    end
                end
            end

            map.external_legs=external_counter-1;
            map.leg_list= leg_list;
            
            map.map = "true";

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