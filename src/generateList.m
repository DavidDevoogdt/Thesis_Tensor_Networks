
%does the same thing but does svd on every O. does not need uniform chain
%also useful for just making the list, untuncated
function mpo_list = generateList(O,M,truncdim,testing)
    if nargin < 4
        testing =0;
    end


    assert( mod(M,2) ==0);
    %assert( M>4);
 
    if iscell(O)
        O_cell = O ;
    else
        O_cell =  truncateO(O,truncdim,testing); %[O_left,O_even,O_odd,O_right]
        
    end
    
    mpo_list = cell(1,M);
 
    O_left =  O_cell{1};
    O_even = O_cell{2};
    O_odd = O_cell{3};
    O_right = O_cell{4};
    
    
    
    %contains O^n
    mpo_list{1} = O_left;
    mpo_list(2:2:M-2) = {O_even};
    mpo_list(3:2:M-1) = {O_odd};
    mpo_list{M} = O_right;


    
end