global d 
d = 2; % d


delta = 1/2; % from xxz model

timestep = 1;

%constuction H tensor final numbering legs: (i1,i2, ..j1,j2)
% basis chosen such that diagonal form persist for diagonalisation
Hij = [delta/4   0           0           0 
       0         -delta/4    1/2         0
       0         1/2         -delta/4    0
       0         0           0           delta/4].*timestep; 
   
H_tensor = reshape(Hij, dimension_vector(d,4));
I_tensor = eye(d);

O_1 = generate_MPO_01( H_tensor,I_tensor,1 );

%some tests


% Make a cell from O. It holds the tensor elements
% every entry is 4d nxdxdxm with n and m the bond dimension for the
% corresponding bond. Dimension x1 at the end not shown by matlab
%this type generates 1--|--1 and 2--|--2 blocks during the expansion
function O = generate_MPO_01( H_tensor,I_tensor, testing)
    global d

    maxIndex = 2;
    O = cell(maxIndex,maxIndex);

    O{1,1} = reshape(  I_tensor, [1,d,d,1] ) ;

    %step 1:
    % 0--|--1--|--0 = exp(H_12) - (0--|--|--0 )  
    N = 1;                  %number of free bonds
    current_max_index = 0;  %used to contract the tensor 

    RHS_Tensor_01 = H_exp(N,H_tensor)- contract_O(N, O ,current_max_index);
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

    RHS_Tensor_11 = H_exp(N,H_tensor) - contract_O(N,O,current_max_index);
    RHS_Matrix_11 = reshape(permute(RHS_Tensor_11, site_ordering_permute(N+1) ),...
                                                  dimension_vector(d,2,[d^2,d^2]) ); 

        
    %todo: this takes the inverse of the O_01 and O_10 matrices, could be done
    % by solving a linear problem or starting from SVD

    O_01_inv = reshape( O{1,2}, [d^2,d^2])^(-1);
    
    O_10_inv = reshape( O{2,1}, [d^2,d^2])^(-1);

    O{2,2} = ncon( { O_01_inv, RHS_Matrix_11, O_10_inv }, { [-1,1], [1,-2,-3,2],[2,-4]}, [1,2]  );
    
    if testing==1
        err = tensor_norm( ncon( { O{1,2}, O{2,2}, O{2,1} }, {[-1,-2,-5,1],[1,-3,-6,2],[2,-4,-7,-8]}, [1,2])...
            -RHS_Tensor_11);
        fprintf("err 11 = %d\n",err);
    end
    
    %step 3:
    % 0--|--1--|--2--|--1--|--0 = exp(H_12+H_23+H_34) - (0--|--|--|--|--0)
    N = 3;                  %number of free bonds
    current_max_index = 1;  %used to contract the tensor 

    RHS_Tensor_12 = H_exp(N,H_tensor) - contract_O(N,O,current_max_index);
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

    RHS_Tensor_22 = H_exp(N,H_tensor) - contract_O(N,O,current_max_index);
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
    
end

%same but does not use 1--|--1 blocks
function O = generate_MPO_02( H_tensor,I_tensor)
    global d

    maxIndex = 2;
    O = cell(maxIndex,maxIndex);

    O{1,1} = reshape(  I_tensor, [1,d,d,1] ) ;
   
    %todo
end


%% all helper fucntions

% return E(H_1_2+..+H_N-1_N) in normal ordering (dim d^N+1, basis first
% upper legs, then lower legs
% so this makes first the tensor T = H x I x I  
%                                  + I x H x I ...
% then reorders, and exponentiates

function H_exp = H_exp(N,H_tensor  )
    global d
    H = zeros( dimension_vector(d, 2*(N+1)) );
    for i = 1:N
        tensor_list = cell(1,N);
        leg_list = cell(1,N);
        
        for j = 1:i-1
            tensor_list{j} = eye(d);
            leg_list{j} = [-j ,-(N+1+j)];
        end
        
        tensor_list{i} = H_tensor;
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

%N number of internal sites
%first upper legs, then lower legs
function T = contract_O(N,O,maxIndex)
    global d
    
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



%helper function to create a 1xn vector  [ left,d,d,..,d,right]
%if left/right are not supplied/0, this is omitted
function p = dimension_vector(d,n,leftright)

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



% this takes a single number and splits it into the composite indices. The
function y = encode_index_array(n,len,d)
    i =1;
    y = zeros(len,1);
    
    while n ~= 0
        [n,r] = readOne(n,d);
        y(i) = r;
        i = i+1;
    end
    
end

% helper function
function [s,r] = readOne(s,d)
    r = rem(s,d);   
    s = (s-r)/d;
end

% this logic was the other way around: define new place a location of old one; 
% function p = site_ordering_permute_old(n)
%         p = zeros(2*n+2,1);
%         p(1)=1;
%         p(2*n+2)=2*n+2;
%         for i = 2:n+1
%             p(i)=2*(i-1);
%         end
%         for i = n+2:2*n+1
%            p(i) = (i-(n+1))*2+1 ;
%         end
% end

% changes from |left i1 i2 ... j1 j2.. right> to |left i1 j1 i2 j2 ...
% right>

function p = site_ordering_permute(n)
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

