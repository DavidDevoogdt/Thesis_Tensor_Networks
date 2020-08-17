
function O_cell_N = mpoProduct(O_cell,N,truncdim,testing,type)
    
    if nargin ==5
       if type==2
           O_cell_N = mpoProduct_assym(O_cell,N,truncdim,testing);
           return
       end
    end

    if nargin < 4
       testing=0;
    end
    
    
    
    O_left = O_cell{1};
    O_even= O_cell{2};
    O_odd = O_cell{3};
    O_right = O_cell{4};
    
    d = size(O_even,2);
    
    oe_dim = size(O_odd,4); %odd even dim
    eo_dim = size(O_odd,1);%even odd dim
    
    %curre t is the running mpo O^N
    O_left_C = O_cell{1};
    O_even_C= O_cell{2};
    O_odd_C = O_cell{3};
    O_right_C = O_cell{4};
    
    for i = 2:N    
        fprintf("adding layer");
        
        oe_dim_C = size(O_odd_C,4); %odd even dim
        eo_dim_C = size(O_odd_C,1);%even odd dim


        Oe_Oo = reshape( ncon( {O_even_C,O_even,O_odd_C,O_odd}, {[-1,-3,4,1],[-2,4,-4,3],[1,-5,2,-7],[3,2,-6,-8]}),...
                                    [ oe_dim*oe_dim_C*d^2, eo_dim*eo_dim_C*d^2 ]);
        [U,S,V] = svds( Oe_Oo, truncdim);
        [a_S_1,td1]=normalisation(S);

        Oe1 = reshape( U*a_S_1, [oe_dim*oe_dim_C,d,d,td1]);
        Oo1 = reshape( S*V'/a_S_1, [td1,d,d,oe_dim*oe_dim_C]);

        Oo1_Oe1 = reshape( ncon( {Oo1,Oe1},{[-1,-2,-3,1] ,[1,-4,-5,-6]}),...
                                        [td1*d^2 ,d^2*td1]);

        [U,S,V] = svds(Oo1_Oe1,truncdim); 
        [a_S_2,td2] = normalisation(S);
        Oo2 = reshape(U*a_S_2, [td1,d,d,td2]);
        Oe2=  reshape(S/a_S_2*V', [td2,d,d,td1]);

        fprintf('tussentijd %s \n', N, datestr(now,'HH:MM:SS.FFF'));

        %fix edges
%         O_left_double = reshape( ncon( {O_left_C,O_left}, {[-1,-3,1,-5],[-2,1,-4,-6]} ),...
%                                     [ 1, d,d, oe_dim*oe_dim_C ]);
% 
%         O_even_1=reshape(Oe1, oe_dim*oe_dim_C,d^2*td1);
%         O_even_2_inv = pinv(  reshape(Oe2, [td2,d^2*td1])  ) ;
%         
%         if testing ==1
%             err= tensor_norm( ncon(  {O_even_1,Oe2 }, {[-1,-2,-4,1],[1,-3,-5,-6]} )-...
%                  ncon(  {O_left_double,Oe1 }, {[-1,-2,-4,1],[1,-3,-5,-6]}  ));
%             fprintf( "beginerr %e \n",err);
%         end
%         
%         begin_O_2 = ncon( {O_left_double, O_even_1,O_even_2_inv}, {[-1,-2,-3,1],[1,2],[2,-4]}  );    
% 
%         if testing ==1
%             err= tensor_norm( ncon(  {begin_O_2,Oe2 }, {[-1,-2,-4,1],[1,-3,-5,-6]} )-...
%                  ncon(  {O_left_double,Oe1 }, {[-1,-2,-4,1],[1,-3,-5,-6]}  ));
%             fprintf( "beginerr %e \n",err);
%         end
% 
%         O_right_double = reshape( ncon( {O_right_C,O_right}, {[-1,-3,1,-5],[-2,1,-4,-6]} ),...
%                                     [ eo_dim*eo_dim_C, d,d,1 ]);
% 
%         O_o=reshape(Oo1, [],eo_dim*eo_dim_C);
%         O_odd_2_inv = pinv(  reshape(Oo2, [td2*d^2,td1])  ) ;
%         end_O_2 = ncon( {O_odd_2_inv,O_o,O_right_double}, {[-1,1],[1,2],[2,-2,-3,-4]}  );    
% 
%         if testing ==1
%             err= tensor_norm( ncon(  {Oo2, end_O_2 }, {[-1,-2,-4,1],[1,-3,-5,-6]} )-...
%                  ncon(  {Oo1, O_right_double }, {[-1,-2,-4,1],[1,-3,-5,-6]}  ));
%             fprintf( "enderr %e \n",err);
%         end


    %%%%%%%%%%%also not accurate
    %search e_2 such that
    % left  --O_odd_C              
    %           |       -- Oe1-- = left_2--Oo2--Oe2--
    % left_c--O_odd    
       
        orig_double = reshape(ncon(  {O_left_C, O_left} , {[-1,-3,1,-5],[-2,1,-4,-6]}),...
                                        [d^2,oe_dim*oe_dim_C ]);
    
        %B
        target_left = reshape( ncon( { orig_double,Oe1}, {[-1,1],[1,-2,-3,-4]}  ),...
                                [d^2,d^2*td1]);
        %A
        target_right = reshape(Oe2 ,[ td1  , d^2* td1]);
        
        %does the same as next command
%         t_left_2 = zeros( d^2, td1);
%         
%         for i = 1:d^2
%             %x*A=B
% 
%             B=target_left(i,:);
%             A=target_right;
%             x = B  /A;
%             t_left_2(i,:) = x;
%             
%         end
        
        %x*A = B
        t_left_2 =  (target_left /target_right) ; %https://nl.mathworks.com/help/matlab/ref/mrdivide.html

        left_2 = reshape(t_left_2, [1,d,d,td1]);
        
        %left_2*target_right-target_left should be about zero
        
        if testing ==1
            err= tensor_norm( ncon( { orig_double,Oe1}, {[-1,1],[1,-2,-3,-4]}  )-...
                 ncon( {t_left_2,Oe2}, {[-1, 1],[1,-2,-3,-4] }  ));
            fprintf( "beginerr %e \n",err);
        end
        
        % same for right
        
        orig_double = reshape(ncon(  {O_right_C, O_right} , {[-1,-3,1,-5],[-2,1,-4,-6]}),...
                                        [eo_dim*eo_dim_C,d^2 ]);
    
        %B
        target_left = reshape( ncon( { Oo1,orig_double}, {[-1,-2,-3,1],[1,-4]}  ),...
                                [d^2*td1,d^2]);
        %A
        target_right = reshape(Oo2 ,[ td1*d^2  ,  td1]);
        
        %x*A = B
        t_right_2 =  (target_right\target_left) ; %https://nl.mathworks.com/help/matlab/ref/mrdivide.html

        right_2 = reshape(t_right_2, [td1,d,d,1]);
        
        %left_2*target_right-target_left should be about zero
        
        if testing ==1
            err= tensor_norm( ncon( { Oo1,orig_double}, {[-1,-2,-3,1],[1,-4]}  )-...
                 ncon( {Oo2,t_right_2}, {[-1,-2,-3, 1],[1,-4] }  ));
            fprintf( "beginerr %e \n",err);
        end
        
       %%%%%%%%%%%%%%%%%%%%%%%% 
        
        O_left_C = left_2;
        O_even_C= Oe2;
        O_odd_C = Oo2;
        O_right_C = right_2;
    
    end
    
    O_cell_N = { O_left_C,O_even_C,O_odd_C,O_right_C};
end




%makes chain of M mpo's next to each other, N layers deep with bond dim =
%truncdim

function mpo_list = mpoProduct_assym(basis_mpo_list,N,truncdim,testing)
    
    if nargin < 6
       testing=0;
    end

   
    mpo_list = basis_mpo_list;
    
    M = size(mpo_list,2);

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
    
    [a_S,td] = normalisation(S);
    
    
    O_left = reshape( U*a_S, [left_bond_dim, d,d, td]);
    O_right = reshape( S/a_S*V', [td,d,d,d_bond_3_right,d_bond_4_right]);
    
end

%s diagonal
function [y,d] = normalisation(S)
    d= min(size(S,1),size(S,2));
    y=0;
    for i = 1:d
       %if S(i,i)>y
       %    y = S(i,i);
       %end
       y=y+S(i,i);
    end
     y=y/d;
end


%just the element wise 2 norm
function norm = tensor_norm(X)
    v = reshape(X,[],1);
    %N=length(v);
    norm = sqrt(  sum(v.^2)  ); 
end