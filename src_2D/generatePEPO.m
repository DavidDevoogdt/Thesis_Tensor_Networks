%for tensors:       ( beta)
%            (alpha)-- O -- (gamma)  = O(i,j,alpha,beta,gamma,delta)
%                    (delta)             1 2   3     4    5     6
% for tensors containig multiple ij's: O numbered from left to right and
% for a given vertical pos from up till down




classdef generatePEPO
    
    properties
        dim
        H_1_tensor
        H_2_tensor
        PEPO
        type
    end
    
    properties (Access = private)
        
    end
    
    methods
        function obj = generatePEPO(d,H_1_tensor,H_2_tensor,max_index)
            obj.dim = d;
            obj.H_1_tensor = H_1_tensor;
            obj.H_2_tensor = H_2_tensor;
            obj.PEPO = cell( max_index+1,max_index+1,max_index+1,max_index+1 );
        end
        
        function H = H_exp(obj,map)
            %map either contains pos_map or (map from parse_map function)   
            map = parse_map(map);

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
                
                H = H + ncon( tensor_list, map.leg_list);
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
                        new_list = [0,0,0,0,0,0,0,0,0,0];
                        
                        new_list( [1,3,5,6,10] )  = index_list_n1([1,2,3,4,6]) ;
                        new_list( [2,4,7,8,9] )  = index_list_n2([1,2,4,5,6]) ;
                        
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
                        
                        %do ij
                        new_list = [0,0,0,0,0,0,0,0,0,0];
                        
                        new_list( [1,3,5,7,8] )  = index_list_n1([1,2,3,4,5]) ;
                        new_list( [2,4,6,9,10] )  = index_list_n2([1,2,3,5,6]) ;
                        
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
        
        function T = contract_network(obj,map, max_index,visualise)
            %generate all index sets for a given configuration
            
            if nargin <4
                visualise = 0;
            end
            
            map = parse_map(map);
            
            if visualise ==1
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
            for n=1: (max_index+1)^map.internal_legs-1
                
                %todo
                tensor_list = cell(1,map.N);

                vect=  encode_index_array(n, map.internal_legs,max_index);
                
                correct_index_set=1;
                
                if visualise==1
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
                        
                        O = obj.PEPO{legs(1)+1,legs(2)+1,legs(3)+1,legs(4)+1};
                        
                        if length(O)==0
                            if visualise==1
                                fprintf("incorrect index set \n");
                            end

                            correct_index_set=0;
                            break;
                        end
                        
                        tensor_list{i} = obj.MPO{ legs };
                    end
                    
                    if ~correct_index_set
                        
                        break;
                    end
                    
                    if visualise==1
                       coor= map.lookup{i} ;
                       new_y = 2+ 2*(coor(1)-1);
                       new_x = 2+ 2*(coor(2)-1);
                       
                       map_copy(new_y,new_x-1) = legs(1);
                       map_copy(new_y,new_x+1) = legs(3);
                       map_copy(new_y-1,new_x) = legs(2);
                       map_copy(new_y+1,new_x) = legs(4);
                       
                    end
                    
                end
                
                if correct_index_set
                    if visualise==1
                        disp(map_copy);
                    end

                    T=T+ncon( tensor_list, map.leg_list);
                end 
            end
        end
    end
end


function map = parse_map(map)
if isfield(map,"pos_map")
    map = create_map(map.pos_map);
else
    if ~isfield(map,"map")
        error("wrong input struct");
    end
end
end

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

map.internal_legs=internal_counter;

%number ij according to number
external_counter=1;

for x =1:n
    for y =1:m
        N1 = pos_map(y,x);
        if N1 ~=0
            leg_list{N1}(1)= -external_counter;
            leg_list{N1}(2)= -(external_counter+N);
            
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