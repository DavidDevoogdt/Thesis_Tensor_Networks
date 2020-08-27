%for tensors:       ( beta)
%            (alpha)-- O -- (gamma)  = O(i,j,alpha,beta,gamma,delta)
%                    (delta)             1 2   3     4    5     6
% for tensors containig multiple ij's: O numbered from left to right and
% for a given vertical pos from up till down


%mxn array with the position 

classdef generatePEPO
    %GENERATEMPO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dim
        H_1_tensor
        H_2_tensor
    end
    
    properties (Access = private)
       
    end
    
    methods
        function obj = generatePEPO(d,H_1_tensor,H_2_tensor)
            %GENERATEMPO Construct an instance of this class
            %   Detailed explanation goes here
            obj.dim = d;
            obj.H_1_tensor = H_1_tensor;
            obj.H_2_tensor = H_2_tensor;
        end
        


        %pos_map
        % array with ones and zeros where contraction mpo are
        %should be conected
        function H = H_exp(obj,pos_map)
            d = obj.dim;
            H_1 = obj.H_1_tensor;
            H_2 = obj.H_2_tensor;
            
            [m,n] = size(pos_map);
            counter = 1;
            for x =1:n
                for y =1:m
                    if pos_map(y,x)==1
                        pos_map(y,x) = counter;
                        counter = counter +1;
                    end
                end
            end

            N = counter-1;

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


            H = zeros( dimension_vector( d,2*N));

            Itensor = reshape(eye(d), [d,d,1,1,1,1]);

            %first do all single site contributions
            for l = 1:N
                tensor_list = cell(1,N);
                tensor_list(:) = {Itensor};
                tensor_list{l} = H_1;

                H = H + ncon( tensor_list, leg_list);
            end

            %do all horizontal H12

           for x =1:n-1
                for y =1:m
                    if pos_map(y,x) ~=0 && pos_map(y,x+1) ~=0
                        n1 = pos_map(y,x);
                        n2 = pos_map(y,x+1);

                        leg_list_copy = leg_list(1,:);

                        index_list_n1 = leg_list{n1};
                        index_list_n2 = leg_list{n2};

                        leg_list_copy(n2) = []; %remove element

                        %do ij
                        new_list = [0,0,0,0,0,0,0,0,0,0];

                        new_list( [1,3,5,6,10] )  = index_list_n1([1,2,3,4,6]) ;
                        new_list( [2,4,7,8,9] )  = index_list_n2([1,2,4,5,6]) ;

                        leg_list_copy{n1} = new_list;

                        tensor_list = cell(1,N-1);
                        tensor_list(:) = {Itensor};
                        tensor_list{n1} = H_2;


                        H=H+ncon( tensor_list,leg_list_copy );
                    end
                end
           end

           %do all vertical H2
           for x =1:n
                for y =1:m-1
                    if pos_map(y,x) ~=0 && pos_map(y+1,x) ~=0
                        n1 = pos_map(y,x);
                        n2 = pos_map(y+1,x);

                        leg_list_copy = leg_list(1,:);

                        index_list_n1 = leg_list{n1};
                        index_list_n2 = leg_list{n2};

                        leg_list_copy(n2) = []; %remove element

                        %do ij
                        new_list = [0,0,0,0,0,0,0,0,0,0];

                        new_list( [1,3,5,7,8] )  = index_list_n1([1,2,3,4,5]) ;
                        new_list( [2,4,6,9,10] )  = index_list_n2([1,2,3,5,6]) ;

                        leg_list_copy{n1} = new_list;

                        tensor_list = cell(1,N-1);
                        tensor_list(:) = {Itensor};
                        tensor_list{n1} = H_2;


                        H=H+ncon( tensor_list,leg_list_copy );
                    end
                end
            end


        end

    end
end
    
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

