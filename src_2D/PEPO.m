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
    end
    
    properties (Access = private)
        
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
            
            O_0000 = expm( obj.H_1_tensor );
            obj.nf = trace(O_0000);
            
            obj.PEPO_cell{1,1,1,1} = reshape(  O_0000/obj.nf , [d,d,1,1,1,1] ) ;

            % 0--|--1--|--0 and all other veriants
            current_max_index = 0;
            
            map = PEPO.create_map([1 1]);
            map_arg = struct( "map",  map);
            
            Tensor_010 = obj.H_exp(map_arg)/(obj.nf^(map.N)) -...
                obj.contract_network(map_arg,current_max_index);
            
            err = eigs(reshape(Tensor_010,[d^(map.N),d^(map.N)]),1);   
            fprintf("010 err %.4e",abs(err));
            
            if abs(err) > 1e-10
                Tensor_010_site = reshape(  permute(Tensor_010, site_ordering_permute(map.N)),...
                                [d^2,d^2] );

                [U,S,V] = svd( Tensor_010_site);


                block_01 = reshape(U, [d,d,d^2]);
                block_10 = permute( reshape(S*V', [d^2,d,d]), [2,3,1]); 

                obj.PEPO_cell{1,1,2,1} =reshape(block_01, [d,d,1,1,d^2,1]);%right
                obj.PEPO_cell{1,1,1,2} =reshape(block_01, [d,d,1,1,1,d^2]);%down

                obj.PEPO_cell{2,1,1,1} =reshape(block_10, [d,d,d^2,1,1,1]);%left
                obj.PEPO_cell{1,2,1,1} =reshape(block_10, [d,d,1,d^2,1,1]);%up


                if obj.testing==1
                    err = ncon( { obj.PEPO_cell{1,1,2,1},  obj.PEPO_cell{2,1,1,1}  },{ [-1,-3,-5,-6,1,-7],[-2,-4,1,-8,-9,-10]  }  )-Tensor_010;
                    fprintf(  "err 010 %.4e \n",  eigs(  reshape(err, [d^map.N,d^map.N]),1) ); 
                end
            end
            
            % create 0--|--1--|--1--|--0 and variants
            current_max_index=1;
            
            map = PEPO.create_map([1 1 1]);
            map_arg = struct( "map",  map);
            
            Tensor_0110 = obj.H_exp(map_arg)/(obj.nf^(map.N))-...
                obj.contract_network(map_arg,current_max_index);
            
            
            err = eigs(reshape(Tensor_0110,[d^(map.N),d^(map.N)]),1);   
            fprintf("0110 err %.4e",abs(err));
            
            if abs(err) > 1e-7
            
                Tensor_0110_site = reshape(  permute(Tensor_0110, site_ordering_permute(map.N)),...
                                [d^2,d^4] );

                left =   reshape(obj.PEPO_cell{1,1,2,1}, [d^2,d^2]) ;
                right = permute( reshape(obj.PEPO_cell{2,1,1,1}, [d^2,d^2]), [2,1]);

                x = left\Tensor_0110_site;
                x = reshape(x, [d^(4),d^(2)]);
                y = x/right;


                block_11 = permute( reshape(y, [d^2,d,d,d^2]), [2,3,1,4]);

                obj.PEPO_cell{2,1,2,1}= reshape(block_11,[d,d,d^2,1,d^2,1]);
                obj.PEPO_cell{2,1,1,2}= reshape(block_11,[d,d,d^2,1,1,d^2]);
                obj.PEPO_cell{1,2,2,1}= reshape(block_11,[d,d,1,d^2,d^2,1]);
                obj.PEPO_cell{1,2,1,2}= reshape(block_11,[d,d,1,d^2,1,d^2]);


                %this is between 2 01 blocks

                map = PEPO.create_map([0 1; 
                                       1 1]);
                map_arg = struct( "map",  map);

                Tensor_0110 = obj.H_exp(map_arg)/(obj.nf^(map.N))-...
                    obj.contract_network(map_arg,current_max_index);

                err = eigs(reshape(Tensor_0110,[d^(map.N),d^(map.N)]),1);   
                fprintf("err %.4e",abs(err));

                Tensor_0110_site = reshape(  permute(Tensor_0110, site_ordering_permute(map.N))....
                                            ,[d^2,d^2,d^2] );

                res= reshape(permute(Tensor_0110_site, [1,3,2]), [d^2,d^4]);

                left =   reshape(obj.PEPO_cell{1,1,2,1}, [d^2,d^2]) ;
                right = permute( reshape(obj.PEPO_cell{1,1,1,2}, [d^2,d^2]), [2,1]);

                x = left\res;
                x = reshape(x, [d^(4),d^(2)]);
                y = x/right;


                block_11 = permute( reshape(y, [d^2,d,d,d^2]), [2,3,1,4]);

                obj.PEPO_cell{2,2,1,1}= reshape(block_11,[d,d,d^2,d^2,1,1]);


                %this is between 2 10 blocks

                map = PEPO.create_map([1 1; 
                                       1 0]);
                map_arg = struct( "map",  map);

                Tensor_0110 = obj.H_exp(map_arg)/(obj.nf^(map.N))-...
                    obj.contract_network(map_arg,current_max_index);

                err = eigs(reshape(Tensor_0110,[d^(map.N),d^(map.N)]),1);   
                fprintf("err %.4e",abs(err));



                Tensor_0110_site = reshape(permute(Tensor_0110, site_ordering_permute(map.N)),...
                                            [d^2,d^2,d^2] );

                res= reshape(permute(Tensor_0110_site, [2,1,3]), [d^2,d^4]);


                left =   reshape(obj.PEPO_cell{1,2,1,1}, [d^2,d^2]) ;
                right = permute( reshape(obj.PEPO_cell{2,1,1,1}, [d^2,d^2]), [2,1]);

                x = left\res;
                x = reshape(x, [d^(4),d^(2)]);
                y = x/right;


                block_11 = permute( reshape(y, [d^2,d,d,d^2]), [2,3,1,4]);

                obj.PEPO_cell{1,1,2,2}= reshape(block_11,[d,d,1,1,d^2,d^2]);

            end
            
            %  correct for one loop
            current_max_index=1;
            map = PEPO.create_map([1 1;
                                   1 1]);
            map_arg = struct( "map",  map);
            
            Tensor_1111 = obj.H_exp(map_arg)/(obj.nf^(map.N))-...
                obj.contract_network(map_arg,current_max_index);
            
            err = eigs(reshape(Tensor_1111,[d^(map.N),d^(map.N)]),1);   
            fprintf("block err %.4e",abs(err));
            
            if abs(err) > 1e-8
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
                        
            end            
        end
        
        function H = H_exp(obj,map)
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
                        n1 = map.num_map(y,x);
                        n2 = map.num_map(y,x+1);
                        
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
                        n1 = map.num_map(y,x);
                        n2 = map.num_map(y+1,x);
                        
                        leg_list_copy = map.leg_list(1,:);
                        
                        index_list_n1 = map.leg_list{n1};
                        index_list_n2 = map.leg_list{n2};
                        
                        leg_list_copy(n2) = []; %remove element
                        
                        for s=1:map.N-1
                            leg_list_copy{s} = leg_list_copy{s}(1:2);
                        end
                        
                        %do ij
                        new_list = [0,0,0,0];
                        
                        new_list( [1,3] )  = index_list_n1([1,2]) ;
                        new_list( [2,4] )  = index_list_n2([1,2]) ;
                        
                        leg_list_copy{n1} = new_list;
                        
                        tensor_list = cell(1,map.N-1);
                        tensor_list(:) = {Itensor};
                        tensor_list{n1} = H_2;
                        
                        
                        H=H+ncon( tensor_list,leg_list_copy );
                    end
                end
            end
            
            H_matrix = expm(reshape(H, [d^(map.N),d^(map.N)]));
            
            H = reshape(H_matrix,dimension_vector(d,2*map.N)  );
        end
        
        function T = contract_network(obj,map, max_index)
            %generate all index sets for a given configuration
                        
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
            T = zeros( dimension_vector( d,2*map.N) );
            
            
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

                    T=T+ncon_optim( tensor_list, map.leg_list);
                end 
            end
        end
        
        function err = calculate_error(obj,map)
            d=obj.dim;
            map = PEPO.parse_map(map);
            map_arg = struct("map",map);
            
            H_matrix=H_exp(obj,map_arg);
            Contraction=obj.contract_network(map_arg,obj.max_index);
            
            a = reshape(  H_matrix, [ d^(map.N), d^(map.N)]);
            b = reshape( Contraction, [ d^(map.N), d^(map.N)]);
            
            trace_a= trace(a);
            
            a = a/trace(a);
            b = b* ( (obj.nf^map.N)/trace_a );
            
            
            err = eigs(a -b,1)  ;
        end
        
    end
    
    methods (Access= private, Static=true)
        function map = parse_map(map)
            if isfield(map,"pos_map")
                map = PEPO.create_map(map.pos_map);
            else
                if isfield(map,"map")
                    map=map.map;
                else
                    error("wrong input struct");
                end
            end
        end
    end
    
    methods (Static=true)
        function map = create_map(pos_map)
            %number the location of operators from up to down and left to
            %right, and create toghether with it a leg_list for ncon

            [m,n] = size(pos_map);
            map.m=m;
            map.n=n;

            map.lookup = {};

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