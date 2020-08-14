classdef generateMPO
    %GENERATEMPO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        d
        H_1_tensor
        H_2_tensor
    end
    
    properties (Access = private) 
        I_tensor
    end
    
    methods
        function obj = generateMPO(d,H_1_tensor,H_2_tensor)
            %GENERATEMPO Construct an instance of this class
            %   Detailed explanation goes here
            obj.d = d;
            obj.H_1_tensor = H_1_tensor;
            obj.H_2_tensor = H_2_tensor;
            obj.I_tensor = eye(d);
        end
        
       
        function MPO = type_01(obj, testing)
            % Make a cell from O. It holds the tensor elements
            % every entry is 4d nxdxdxm with n and m the bond dimension for the
            % corresponding bond. dimension x1 at the end not shown by matlab
            %this type generates 1--|--1 and 2--|--2 blocks during the expansion
            d = obj.d;

            maxIndex = 2;
            O = cell(maxIndex,maxIndex);

            O{1,1} = reshape(  obj.I_tensor, [1,d,d,1] ) ;

            %step 1:
            % 0--|--1--|--0 = exp(H_12) - (0--|--|--0 )  
            N = 1;                  %number of free bonds
            current_max_index = 0;  %used to contract the tensor 

            RHS_Tensor_01 = H_exp(obj,N)- contract_O(N, O ,current_max_index,d);
            RHS_Matrix_01 = reshape( permute( RHS_Tensor_01 , site_ordering_permute(N+1) ),...
                                                        [d^2,d^2] ); %ready to svd

            [U,S,V] = svd(RHS_Matrix_01);
            sqrt_S = sqrt(S); %for symmetric split, not really necesary

            O{0+1,1+1} = reshape( U*sqrt_S, [1,d,d,d^2]);
            O{1+1,0+1} = reshape( sqrt_S* V', [d^2,d,d,1]);

            if testing==1
                err = tensor_norm( ncon( {O{0+1,1+1},O{1+1,0+1}}, {[-1,-2,-4,1],[1,-3,-5,-6]}, [1])-RHS_Tensor_01);
                fprintf("err 01 = %d\n",err);
            end

            %step 2 :
            % 0--|--1--|--1--|--0 = exp(H12+H23)- (0--|--|--|--0) 
            N = 2;                  %number of free bonds
            current_max_index = 1;  %used to contract the tensor 

            RHS_Tensor_11 =H_exp(obj,N) - contract_O(N,O,current_max_index,d);
            RHS_Matrix_11 = reshape(permute(RHS_Tensor_11, site_ordering_permute(N+1) ),...
                                                          [d^2,d,d,d^2] ); 


            %todo: this takes the inverse of the O_01 and O_10 matrices, could be done
            % by solving a linear problem or starting from SVD

            O_01_inv = reshape( O{1,2}, [d^2,d^2])^(-1);
            O_10_inv = reshape( O{2,1}, [d^2,d^2])^(-1);
            O{2,2} = ncon( { O_01_inv, RHS_Matrix_11, O_10_inv }, { [-1,1], [1,-2,-3,2],[2,-4]}, [1,2]  ); 
%             
%             O_01 = reshape( O{1,2}, [d^2,d^2]);
% 
%             O_10 = reshape( O{2,1}, [d^2,d^2]);
%             
%             temp =  O_01\RHS_Matrix_11;  %https://nl.mathworks.com/help/matlab/ref/mldivide.html
%             O_11 = temp/O_10;
            

            %O{2,2} = reshape( O_11, [d^2,d,d,d^2]);

            if testing==1
                err = tensor_norm( ncon( { O{1,2}, O{2,2}, O{2,1} }, {[-1,-2,-5,1],[1,-3,-6,2],[2,-4,-7,-8]}, [1,2])...
                    -RHS_Tensor_11);
                fprintf("err 11 = %d\n",err);
            end

            %step 3:
            % 0--|--1--|--2--|--1--|--0 = exp(H_12+H_23+H_34) - (0--|--|--|--|--0)
            N = 3;                  %number of free bonds
            current_max_index = 1;  %used to contract the tensor 

            RHS_Tensor_12 = H_exp(obj,N)- contract_O(N,O,current_max_index,d);
            RHS_Matrix_12 = reshape( permute(RHS_Tensor_12, site_ordering_permute(N+1) ),...
                                                   dimension_vector(d^2,4));%group per ij index


            O_12_21 = ncon( { O_01_inv, RHS_Matrix_12, O_10_inv }, { [-1,1], [1,-2,-3,2],[2,-4]}, [1,2]  );
            O_12_21_svd = reshape(O_12_21, [d^4,d^4]);

            [U,S,V] = svd(O_12_21_svd);
            sqrt_S = sqrt(S); %for symmetric split, not really necesary

            O{1+1,2+1} = reshape( U*sqrt_S, [d^2,d,d,d^4]);
            O{2+1,1+1} = reshape( sqrt_S* V', [d^4,d,d,d^2]);

            if testing==1
                err = tensor_norm( ncon( {O{1,2},O{2,3},O{3,2},O{2,1}}, {[-1,-2,-6,1],[1,-3,-7,2],[2,-4,-8,3],[3,-5,-9,-10]}, [1,2,3])...
                    -RHS_Tensor_12);
                fprintf("err 12 = %d\n",err);
            end

            %step 4:
            % 0--|--1--|--2--|--2--|--1--|--0 = exp(H_12+H_23+H_34+H_45) - (0--|--|--|--|--|--0)
            N = 4;                  %number of free bonds
            current_max_index = 2;  %used to contract the tensor 

            RHS_Tensor_22 = H_exp(obj,N) - contract_O(N,O,current_max_index,d);
            RHS_Tensor_22_site = reshape( permute(RHS_Tensor_22, site_ordering_permute(N+1)),...
                                    dimension_vector(d^2,5)); %group per ij index

            % first apply inverses of O_01 and O_02
            RHS_Tensor_22_site_stripped = ncon( {O_01_inv, RHS_Tensor_22_site, O_10_inv}, { [-1,1],[1,-2,-3,-4,2],[2,-5]},[1,2]);


            O_12_inv = (reshape( O{1+1,2+1}, [d^4,d^4]))^-1;  % (O_12, (alpha,i,j) beta)^-1
            O_21_inv = (reshape( O{2+1,1+1}, [d^4,d^4]))^-1;

            RHS_Tensor_22_site_stripped_reshaped = reshape( RHS_Tensor_22_site_stripped, [d^4,d,d,d^4]);

            O{2+1,2+1} = ncon(  {O_12_inv,RHS_Tensor_22_site_stripped_reshaped,O_21_inv},{ [-1,1],[1,-2,-3,2],[2,-4]},[1,2]);


            if testing==1
                err = tensor_norm( ncon( {O{1,2},O{2,3},O{3,3},O{3,2},O{2,1}}, {[-1,-2,-7,1],[1,-3,-8,2],[2,-4,-9,3],[3,-5,-10,4],[4,-6,-11,-12]}, [1,2,3,4])...
                    -RHS_Tensor_22);
                fprintf("err 12 = %d\n",err);
            end

            %reshape such that the ruslting mpo is one matrix instead of several 
            %virtual levels
            %each seperate element is in ordering  upper legs.. lower legs
            

            MPO = mpo_cell_2_matrix(O,maxIndex,d);
            
            %seems to work
            if testing==1
                left =zeros(1,21);
                left(1) = 1;
                right =zeros(21,1);
                right(1) = 1;
                
                err1 = tensor_norm( ncon( {left,MPO,right},{[-1,1],[1,-2,-3,2],[2,-4]})-contract_O(0,O,2,d));
                err2 = tensor_norm( ncon( {left,MPO,MPO,right},{[-1,1],[1,-2,-4,2],[2,-3,-5,3],[3,-6]})-contract_O(1,O,2,d));
                err3 = tensor_norm( ncon( {left,MPO,MPO,MPO,right},{[-1,1],[1,-2,-5,2],[2,-3,-6,3],[3,-4,-7,4],[4,-8]})-contract_O(2,O,2,d));
                fprintf("matrixconversion  err0 %d err1 %d err3 %d \n",err1,err2,err3);
            end
            

            
         
            
        end
        
        
        function MPO = type_02(obj, testing)
            % Make a cell from O. It holds the tensor elements
            % every entry is 4d nxdxdxm with n and m the bond dimension for the
            % corresponding bond. dimension x1 at the end not shown by matlab
            %this type generates no 1--|--1 and 2--|--2 blocks
            d = obj.d;

            maxIndex = 2;
            O = cell(maxIndex,maxIndex);

            O{1,1} = reshape(  obj.I_tensor, [1,d,d,1] ) ;

            %step 1:
            % 0--|--1--|--0 = exp(H_12) - (0--|--|--0 )  
            N = 1;                  %number of free bonds
            current_max_index = 0;  %used to contract the tensor 

            RHS_Tensor_01 = H_exp(obj,N)- contract_O(N, O ,current_max_index,d);
            RHS_Matrix_01 = reshape( permute( RHS_Tensor_01 , site_ordering_permute(N+1) ),...
                                                        [d^2,d^2] ); %ready to svd

            [U,S,V] = svd(RHS_Matrix_01);
            sqrt_S = sqrt(S); %for symmetric split, not really necesary

            O{0+1,1+1} = reshape( U*sqrt_S, [1,d,d,d^2]);
            O{1+1,0+1} = reshape( sqrt_S* V', [d^2,d,d,1]);

            if testing==1
                err = tensor_norm( ncon( {O{0+1,1+1},O{1+1,0+1}}, {[-1,-2,-4,1],[1,-3,-5,-6]}, [1])-RHS_Tensor_01);
                fprintf("err 01 = %d\n",err);
            end

            %step 2 :
            % 0--|--1--|--2--|--0 = exp(H12+H23)- (0--|--|--|--0) 
            N = 2;                  %number of free bonds
            current_max_index = 1;  %used to contract the tensor 

            RHS_Tensor_12 = H_exp(obj,N) - contract_O(N,O,current_max_index,d);
            RHS_Matrix_12 = reshape(permute(RHS_Tensor_12, site_ordering_permute(N+1) ),...
                                                          [d^2,d^2,d^2] ); 


            %todo: this takes the inverse of the O_01 and O_10 matrices, could be done
            % by solving a linear problem or starting from SVD

            O_01_inv = reshape( O{1,2}, [d^2,d^2])^(-1);
            
            RHS_Tensor_120 = reshape( ncon( { O_01_inv, RHS_Matrix_12 }, { [-1,1], [1,-2,-3,-4]}  ),...
                                       [d^4,d^2] );
                
            [U,S,V] = svd(RHS_Tensor_120);
            %sqrt_S = sqrt(S); %for symmetric split, not really necesary

            O{1+1,2+1} = reshape( U, [d^2,d,d,d^4]);
            O{2+1,0+1} = reshape( S*V', [d^4,d,d,1]);
                                   
            if testing==1
                err = tensor_norm( ncon( { O{1,2}, O{2,3}, O{3,1} }, {[-1,-2,-5,1],[1,-3,-6,2],[2,-4,-7,-8]}, [1,2])...
                    -RHS_Tensor_12);
                fprintf("err 0120 = %d\n",err);
            end
           
            %step 2 :
            % 0--|--1--|--2--|--1--|--0 = exp(H12+H23+h34)- (0--|--|--|--|--0) 
            N = 3;                  %number of free bonds
            current_max_index = 1;  %used to contract the tensor 

            RHS_Tensor_121 = H_exp(obj,N) - contract_O(N,O,current_max_index,d);
            RHS_Matrix_121 = reshape(permute(RHS_Tensor_121, site_ordering_permute(N+1) ),...
                                                          [d^2,d^2,d^2,d^2] ); 


            %todo: this takes the inverse of the O_01 and O_10 matrices, could be done
            % by solving a linear problem or starting from SVD

            
            O_12_inv = reshape( O{2,3}, [d^4,d^4])^(-1);
            O_10_inv = reshape( O{2,1}, [d^2,d^2])^(-1);
            
            temp = reshape(ncon( {O_01_inv, RHS_Matrix_121, O_10_inv }, { [-1,2] [2,-2,-3,3],[3,-4]}  ) ,...
                                       [d^4,d,d,d^2 ] );
                                   
            O{3,2} =  ncon( {O_12_inv,temp}, {[-1,1],[1,-2,-3,-4]}  ) ;                    
          
                                   
            if testing==1
                err = tensor_norm( ncon( { O{1,2}, O{2,3}, O{3,2},O{2,1} }, {[-1,-2,-6,1],[1,-3,-7,2],[2,-4,-8,3] ,[3,-5,-9 ,-10 ]} )...
                    -RHS_Tensor_121);
                fprintf("err 01210 = %d\n",err);
            end
            
            
            
            
            MPO = mpo_cell_2_matrix(O,maxIndex,d);
            
            %seems to work
            if testing==1
                left =zeros(1,21);
                left(1) = 1;
                right =zeros(21,1);
                right(1) = 1;
                
                err1 = tensor_norm( ncon( {left,MPO,right},{[-1,1],[1,-2,-3,2],[2,-4]})-contract_O(0,O,2,d));
                err2 = tensor_norm( ncon( {left,MPO,MPO,right},{[-1,1],[1,-2,-4,2],[2,-3,-5,3],[3,-6]})-contract_O(1,O,2,d));
                err3 = tensor_norm( ncon( {left,MPO,MPO,MPO,right},{[-1,1],[1,-2,-5,2],[2,-3,-6,3],[3,-4,-7,4],[4,-8]})-contract_O(2,O,2,d));
                fprintf("matrixconversion  err0 %d err1 %d err3 %d \n",err1,err2,err3);
            end
            
            
        end
        

    end
    
    methods (Access=private)
        

        function H_exp = H_exp(obj,N)
            % return E(H_1_2+..+H_N-1_N) in normal ordering (dimension d^N+1, basis first
            % upper legs, then lower legs
            % so this makes first the tensor T = H     x I x I  
            %                                  + I x H     x I ...
            %
            %                                  + S x I x I x I       
            %                                  + I x S x I x I ... 
            % with S the single site operator and H the 2 site one
            % then reorders, and exponentiates
            d = obj.d;
            
            H = zeros( dimension_vector(d, 2*(N+1)) );
            
            %1 site operator
            for i = 1:N+1
                tensor_list = cell(1,N+1);
                leg_list = cell(1,N+1);

                for j = 1:i-1
                    tensor_list{j} = eye(d);
                    leg_list{j} = [-j ,-(N+1+j)];
                end

                tensor_list{i} = obj.H_1_tensor;
                leg_list{i} = [-i, -(N+1+i)  ];

                for j = i+1:N+1
                    tensor_list{j} = eye(d);
                    leg_list{j} = [-(j) ,-(N+1+j)];
                end

                H_i = ncon( tensor_list, leg_list, []);
                H = H + H_i;
            end
            
            
            %2 site operator
            for i = 1:N
                tensor_list = cell(1,N);
                leg_list = cell(1,N);

                for j = 1:i-1
                    tensor_list{j} = eye(d);
                    leg_list{j} = [-j ,-(N+1+j)];
                end

                tensor_list{i} = obj.H_2_tensor;
                 leg_list{i} = [-i, -(i+1), -(N+1+i), -(N+i+2)  ];

                for j = i+1:N
                    tensor_list{j} = eye(d);
                    leg_list{j} = [-(j+1) ,-(N+j+2)];
                end

                H_i = ncon( tensor_list, leg_list, []);
                H = H + H_i;
            end

            H_exp = reshape(  expm( reshape(H, [d^(N+1),d^(N+1)])), dimension_vector(d,2*(N+1),[1,1] ));


        end
    end

end



function T = contract_O(N,O,maxIndex,d)
    %N number of internal sites
    %returns first upper legs, then lower legs
    
    %global d

    T = zeros( dimension_vector(d,2*(N+1),[1,1]));
    %generate all combinations of internam indices
    for i= 0: (maxIndex+1)^N-1
        full_vect = [0;encode_index_array(i,N,maxIndex+1);0];

        correct_index_set = 1;

        O_tensors =  cell(1,N+1);
        %Take correct O's for contraction with ncon
        for j = 1:N+1 
            O_j = O{  full_vect(j)+1 ,full_vect(j+1)+1};
            if length(O_j)==0
                correct_index_set = 0;
                break;
            end
            O_tensors{j} = O_j;
        end


        if correct_index_set == 1
            legs = cell(1,N+1);
            for t=1:N 
                legs{t} = zeros(1,4);
            end

            legs{1}(1) = -1;%first index
            legs{N+1}(4) = -(  2*(N+1) +2); % last one

            for t=1:N %assign all indices to sum over
                legs{t}(4) = t;
                legs{t+1}(1) = t;
            end

            for t=1:N+1 %assign final places for i_n and j_n
                legs{t}(2) = -(t+1); 
                legs{t}(3) = -(N+1 +t+1);
            end

            seq = zeros(1,N);
            for t=1:N %todo think harder about contraction order
                seq(t)=t;
            end

            T_j = ncon (O_tensors,legs, seq);
            T = T + T_j;
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




function y = encode_index_array(n,len,d)
    % this takes a single number and splits it into the composite indices. The
    i =1;
    y = zeros(len,1);

    while n ~= 0
        [n,r] = readOne(n,d);
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
    p = zeros(2*n+2,1);
    p(1)=1;
    p(2*n+2)=2*n+2;
    for i = 1:n
        p(2*i)=i+1;
    end
    for i = 1:n
       p(2*i+1) = n+i+1 ;
    end
end

%just the element wise 2 norm
function norm = tensor_norm(X)
    v = reshape(X,[],1);
    %N=length(v);
    norm = sqrt(  sum(v.^2)  ); 
end


function T = mpo_cell_2_matrix(O,maxIndex,d)


    totaldimension = geomSum(d^2,maxIndex-1);

    T = zeros(totaldimension,d,d,totaldimension);

    %reassemble into a large tensor without cells
    for i = 0:maxIndex
        for j = 0:maxIndex
            start_i= geomSum(d^2,i-1);
            end_i =  start_i+d^(2*i)-1;
            start_j= geomSum(d^2,j-1);
            end_j =  start_j+d^(2*j)-1;

            O_ij = O{i+1,j+1};
            if length(O_ij) ~= 0
                T( start_i+1:end_i+1 ,:,:, start_j+1:end_j+1 ) = O{i+1,j+1};
            end
        end
    end

    function z =geomSum(x,n)
       z= (x^(n+1)-1)/(x-1);
    end


end

