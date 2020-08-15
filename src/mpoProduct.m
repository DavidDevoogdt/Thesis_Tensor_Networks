%makes chain of M mpo's next to each other, N layers deep with bond dim =
%truncdim

function mpo_list = mpoProduct(O,M,N,truncdim,testing)
    if nargin < 5
       testing=0;
    end

    d = size(O,2);
    chi = size(O,1);
    
    end_vector = zeros(1,chi);
    end_vector(1)=1;
    
    mpo_left = ncon(  {end_vector,O}, {[-1,1],[1,-2,-3,-4]});
    mpo_right = ncon(  {O,end_vector}, {[-1,-2,-3,1],[-4,1]});
    
    mpo_list = cell(1,M);
    
    
    %contains O^n
    mpo_list(2:M-1) = {O};
    mpo_list{1} = mpo_left;
    mpo_list{M} = mpo_right;

    
    basis_mpo_list = mpo_list; %list of O to apply
    
    if chi > truncdim %pretruncate
        right = mpo_list{1};
        for i = 1:M-1
            
            con = ncon( {right,mpo_list{i+1}},{  [-1,-2,-3,1],[1,-4,-5,-6] });
            size_vect =size(con);
            con = reshape(con,  size_vect(1)*size_vect(2)*size_vect(3), []);

            [U,S,V] = svds( con,truncdim);
            %from this paper: scale S and U https://arxiv.org/pdf/1611.02498.pdf
            [a_S,td] = average(S);
          
            left = reshape( U*a_S,  size_vect(1) , d,d, td);
            right = reshape( S/a_S*V', td,d,d, []);
            
            basis_mpo_list{i} = left;
        end
        basis_mpo_list{M} = right;
        
        mpo_list = basis_mpo_list;
    end
    
   
    
    for n = 2:N
        
        [left, right] = mpo_product_4(mpo_list{1}, basis_mpo_list{1} ,mpo_list{2},basis_mpo_list{2},truncdim);
        mpo_list{1} = left;
        
        for i = 2:M-2
            [left,right] = mpo_product_3(right, mpo_list{i+1} ,basis_mpo_list{i+1},truncdim);
            mpo_list{i} = left;
        end

        [left,right] = mpo_product_3(right,mpo_list{M},basis_mpo_list{M},truncdim);
        mpo_list{M-1} = left;
        mpo_list{M} = right;
        if testing==1
           fprintf("added layer %d\n",n); 
        end
    end
    
end

%connects the 4 mpo's to 2 mpo's with bond dim = truncdim  = td
%    d|      |
% -- O_1 -- O_3 --                  |  td  |
%     |      |          =>      -- O_5 -- O_6 --
% -- O_2 -- O_4 --                  |      |
%     |      |   
function [O_left,O_right] = mpo_product_4(mpo_1,mpo_2,mpo_3,mpo_4,truncdim)

    con =  ncon( {mpo_1,mpo_2}, { [-1,-3,1,-5],[-2,1,-4,-6]});
    size_vect = size(con);
    mpo_A = reshape(con, [size_vect(1)*size_vect(2), size_vect(3),size_vect(4),size_vect(5),size_vect(6)]);

    [O_left,O_right] = mpo_product_3(mpo_A,mpo_3,mpo_4,truncdim);

    
end

%connects the 4 mpo's to 2 mpo's with bond dim = truncdim  = td
%    d|      |
%        -- O_3 --                  |  td  |  --
% -- O_1     |          =>      -- O_5 -- O_6 
%        -- O_4 --                  |      |  --
%     |      |   

function [O_left,O_right] = mpo_product_3(mpo_1,mpo_3,mpo_4,truncdim)
 
    

    d = size(mpo_1,2); %physical dimension
   
    left_bond_dim = size(mpo_1,1); 
    
    d_bond_3_right = size(mpo_3,4);
    d_bond_4_right = size(mpo_4,4);
    
    right_bond_dim = d_bond_3_right*d_bond_4_right;
    
    assert(size(mpo_1,4)==size(mpo_3,1));
    assert(size(mpo_1,5)==size(mpo_4,1));
    
    
    sum = ncon( {mpo_1,mpo_3,mpo_4}, { [-1,-2,-3,1,2],[1,-4,3,-6],[2,3,-5,-7]});
     
    [U,S,V] = svds(  reshape( sum, [left_bond_dim* d^2,d^2*right_bond_dim]),truncdim);
    
    [a_S,td] = average(S);
    
    
    O_left = reshape( U*a_S, [left_bond_dim, d,d, td]);
    O_right = reshape( S/a_S*V', [td,d,d,d_bond_3_right,d_bond_4_right]);
    
end

%s diagonal
function [y,d] = average(S)
    d= min(size(S,1),size(S,2));
    y=0;
    for i = 1:d
       y=y+S(i,i) ;
    end

end
