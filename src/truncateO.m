%generates list O-O-O-O for given O, with max bond dim truncdim
% D--O--D--O-- -> D---O_left--d--O_right--D
% d--O_right--D---O_left--d -> 
% O_cell = O_left,O_even,O_odd,O_right}
function O_cell = truncateO(O,truncdim,testing)
    chi = size(O,1);
    d =size(O,2);

%for future impl
    end_vector = zeros(1,chi);
    end_vector(1) = 1;

    left = end_vector;
    right = left';
    
   
    
    O_odd = O;
    O_even = O;

    if nargin < 3
        testing=0;
    end


    if truncdim==-1
         O_left_2= ncon( {left,O}, {  [-1,1],[1,-2,-3,-4] }  );
         O_right_2 = ncon( {O,right}, {[-1,-2,-3,1],[1,-4]}  );
         O_cell = {O_left_2,O_odd,O_even,O_right_2};
         return
    end

    Oe_Oo = ncon(  {O_even,O_odd},{[-1,-2,-3,1],[1,-4,-5,-6]} );
    [U,S_1,V] = svds(  reshape( Oe_Oo, [chi*d^2,d^2*chi ]),truncdim );

    [a_S_1,td_1] = normalisation(S_1);

    O_even_1 = reshape( U*a_S_1, [chi, d,d,td_1]);
    O_odd_1 = reshape( S_1/a_S_1*V', [td_1,d,d,chi]);


    [U2,S_2,V2] = svds(reshape(ncon( {O_odd_1,O_even_1}, {[-1,-2,-3,1],[1,-4,-5,-6]}), [ td_1*d^2, d^2*td_1  ]  ),truncdim );

    [a_S_2,td_2] = normalisation(S_2);

    O_odd_2 = reshape( U2*a_S_2, [td_1, d,d,td_2]);
    O_even_2 = reshape( S_2/a_S_2*V2', [td_2,d,d,td_1]);

    % recalculate the left and right endvectors 
    
    %good enough for here
    %search e_2 such that
    % left--O_odd--O_even_1-- = left_2--O_odd_2--O_even_2--

    target_left = reshape( ncon( {left,O_odd,O_even_1}, {[-1,1],[1,-2,-3,2],[2,-4,-5,-6] }  ),...
                            [1,d^4*td_1]);
    
    target_right = reshape( ncon( {O_odd_2,O_even_2}, {[-1,-2,-3,2],[2,-4,-5,-6] }  ),...
                            [ td_1  , d^4*td_1]);
    
    left_2 =  (target_left /target_right) ; %https://nl.mathworks.com/help/matlab/ref/mrdivide.html
    

    if testing ==1
        err= tensor_norm(ncon( {left,O_odd,O_even_1}, {[-1,1],[1,-2,-3,2],[2,-4,-5,-6] }  )-...
             ncon( {left_2,O_odd_2,O_even_2}, {[-1,1],[1,-2,-3,2],[2,-4,-5,-6] }  ));
        fprintf( "beginerr %e \n",err);
    end
    
    
   
    target_left = reshape( ncon( {O_odd_1,O_even,right}, {[-1,-2,-3,2],[2,-4,-5,1],[1,-6] }  ),...
                            [d^4*td_1,1]);
    
    target_right = reshape( ncon( {O_odd_2,O_even_2}, {[-1,-2,-3,2],[2,-4,-5,-6] }  ),...
                            [ td_1*d^4  , td_1]);
    
   
    right_2 =  ( target_right \target_left) ; %https://nl.mathworks.com/help/matlab/ref/mrdivide.html

    if testing ==1
        err= tensor_norm(ncon( {O_odd_1,O_even,right}, {[-1,-2,-3,2],[2,-4,-5,1],[1,-6] }  )-...
             ncon( {O_odd_2,O_even_2,right_2}, {[-1,-2,-3,2],[2,-4,-5,1],[1,-6] }  ));
        fprintf( "beginerr %e \n",err);
    end

    
    left = ncon( {left_2,O_odd_2}, {  [-1,1],[1,-2,-3,-4] }  );
    
    right = ncon( {O_even_2,right_2}, {[-1,-2,-3,1],[1,-4]}  );

    %assign return values
    O_cell = {left ,O_even_2,O_odd_2,right};
end

    
%s diagonal
function [y,d] = normalisation(S)
    d= min(size(S,1),size(S,2));
    y=0;
    for i = 1:d
       if S(i,i)>y
           y = S(i,i);
       end
    end
     y=y/d;
end

%just the element wise 2 norm
function norm = tensor_norm(X)
    v = reshape(X,[],1);
    %N=length(v);
    norm = sqrt(  sum(v.^2)  ); 
end