classdef generateMPO
    %GENERATEMPO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dim
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
            obj.dim = d;
            obj.H_1_tensor = H_1_tensor;
            obj.H_2_tensor = H_2_tensor;
            obj.I_tensor = eye(d);
        end
        
        
        function [normalisation_factor,MPO] = type_01(obj, testing)
            % Make a cell from O. It holds the tensor elements
            % every entry is 4d nxdxdxm with n and m the bond dimension for the
            % corresponding bond. dimension x1 at the end not shown by matlab
            %this type generates 1--|--1 and 2--|--2 blocks during the expansion
            d = obj.dim;
            
            maxIndex = 2;
            O = cell(maxIndex+1,maxIndex+1);
            
            O_11_unnormalised = expm( obj.H_1_tensor );
            normalisation_factor = trace(O_11_unnormalised);
            
            O{1,1} = reshape(  O_11_unnormalised/normalisation_factor , [1,d,d,1] ) ;
            
            
            
            %step 1:
            % 0--|--1--|--0 = exp(H_12) - (0--|--|--0 )
            N = 1;                  %number of free bonds
            current_max_index = 0;  %used to contract the tensor
            
            RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1)) - contract_O(N, O ,current_max_index,d);
            RHS_Matrix = reshape( permute( RHS_Tensor , site_ordering_permute(N+1) ),...
                [d^2,d^2] ); %ready to svd
            
            [U,S,V] = svd(RHS_Matrix);
            %sqrt_S = sqrt(S); %for symmetric split, not really necesary
            
            a_S = average(S);
            
            O{0+1,1+1} = reshape( U*a_S, [1,d,d,d^2]);
            O{1+1,0+1} = reshape( S/a_S* V', [d^2,d,d,1]);
            
            if testing==1
                err = tensor_norm( ncon( {O{0+1,1+1},O{1+1,0+1}}, {[-1,-2,-4,1],[1,-3,-5,-6]}, [1])-RHS_Tensor);
                fprintf("err 010 = %d\n",err);
            end
            
            %step 2 :
            % 0--|--1--|--1--|--0 = exp(H12+H23)- (0--|--|--|--0)
            N = 2;                  %number of free bonds
            current_max_index = 1;  %used to contract the tensor
            
            RHS_Tensor =H_exp(obj,N)/(normalisation_factor^(N+1)) - contract_O(N,O,current_max_index,d);
            RHS_Matrix = reshape(permute(RHS_Tensor, site_ordering_permute(N+1) ),...
                [d^2,d,d,d^2] );
            
            
            %search x st
            % left*x = res
            % y*right = x
            left = reshape( O{1,2}, [d^2,d^2]);
            right = reshape( O{2,1}, [d^2,d^2]);
            res = reshape( RHS_Matrix, [d^2,d^4]);
            
            x = left\res;
            x = reshape(x, [d^4,d^2]);
            y = x/right;
            
            
            O{2,2} = reshape(y, [d^2,d,d,d^2 ]);
            
            
            if testing==1
                err = tensor_norm( ncon( { O{1,2}, O{2,2}, O{2,1} }, {[-1,-2,-5,1],[1,-3,-6,2],[2,-4,-7,-8]}, [1,2])...
                    -RHS_Tensor);
                fprintf("err 0110 = %d\n",err);
            end
            
            %step 3:
            % 0--|--1--|--2--|--1--|--0 = exp(H_12+H_23+H_34) - (0--|--|--|--|--0)
            N = 3;                  %number of free bonds
            current_max_index = 1;  %used to contract the tensor
            
            RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1))- contract_O(N,O,current_max_index,d);
            RHS_Matrix = reshape( permute(RHS_Tensor, site_ordering_permute(N+1) ),...
                dimension_vector(d^2,4));%group per ij index
            
            %search x st
            % left*x = res
            % y*right = x
            left = reshape( O{1,2}, [d^2,d^2]);
            right = reshape( O{2,1}, [d^2,d^2]);
            res = reshape( RHS_Matrix, [d^2,d^6]);
            
            x = left\res;
            x = reshape(x, [d^6,d^2]);
            y = x/right;
            
            new_parts = reshape(y, [d^4,d^4]);
            
            [U,S,V] = svd(new_parts);
            a_S = average(S);
            
            
            O{1+1,2+1} = reshape( U*a_S, [d^2,d,d,d^4]);
            O{2+1,1+1} = reshape( S/a_S* V', [d^4,d,d,d^2]);
            
            if testing==1
                err = tensor_norm( ncon( {O{1,2},O{2,3},O{3,2},O{2,1}}, {[-1,-2,-6,1],[1,-3,-7,2],[2,-4,-8,3],[3,-5,-9,-10]}, [1,2,3])...
                    -RHS_Tensor);
                fprintf("err 01210 = %d\n",err);
            end
            
            %step 4:
            % 0--|--1--|--2--|--2--|--1--|--0 = exp(H_12+H_23+H_34+H_45) - (0--|--|--|--|--|--0)
            N = 4;                  %number of free bonds
            current_max_index = 2;  %used to contract the tensor
            
            RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1)) - contract_O(N,O,current_max_index,d);
            RHS_Matrix = reshape( permute(RHS_Tensor, site_ordering_permute(N+1)),...
                dimension_vector(d^2,5)); %group per ij index
            
            %search x st
            % left*x = res
            % y*right = x
            left = reshape( ncon( {O{1,2},O{2,3}},{[-1,-2,-3,1],[1,-4,-5,-6]}),...
                [d^4,d^4]);
            right = reshape( ncon( {O{3,2},O{2,1} },{[-1,-2,-3,1],[1,-4,-5,-6]}),...
                [d^4,d^4]);
            res = reshape( RHS_Matrix, [d^4,d^6]);
            
            x = left\res;
            x = reshape(x, [d^6,d^4]);
            y = x/right;
            
            
            O{2+1,2+1} = reshape(y, [d^4,d,d,d^4]);
            
            
            if testing==1
                err = tensor_norm( ncon( {O{1,2},O{2,3},O{3,3},O{3,2},O{2,1}}, {[-1,-2,-7,1],[1,-3,-8,2],[2,-4,-9,3],[3,-5,-10,4],[4,-6,-11,-12]}, [1,2,3,4])...
                    -RHS_Tensor);
                fprintf("err 012210 = %d\n",err);
            end
            
            %reshape such that the ruslting mpo is one matrix instead of several
            %virtual levels
            %each seperate element is in ordering  upper legs.. lower legs
            
            
            MPO = mpo_cell_2_matrix(O,maxIndex,d);
            
            %seems to work
            %             if testing==1
            %                 left =zeros(1,21);
            %                 left(1) = 1;
            %                 right =zeros(21,1);
            %                 right(1) = 1;
            %
            %                 err1 = tensor_norm( ncon( {left,MPO,right},{[-1,1],[1,-2,-3,2],[2,-4]})-contract_O(0,O,2,d));
            %                 err2 = tensor_norm( ncon( {left,MPO,MPO,right},{[-1,1],[1,-2,-4,2],[2,-3,-5,3],[3,-6]})-contract_O(1,O,2,d));
            %                 err3 = tensor_norm( ncon( {left,MPO,MPO,MPO,right},{[-1,1],[1,-2,-5,2],[2,-3,-6,3],[3,-4,-7,4],[4,-8]})-contract_O(2,O,2,d));
            %                 fprintf("matrixconversion  err0 %d err1 %d err3 %d \n",err1,err2,err3);
            %             end
            %
            
            
            
            
        end
        
        function [normalisation_factor,MPO] = type_04(obj,order,opts)
            % this type generates n--|--n blocks during the expansion
            %order n means explicit calculation op to order free bonds
            d = obj.dim;
            
            
            p = inputParser;
            addParameter(p,'method',"svd")
            addParameter(p,'testing',0)
            parse(p,opts)
            
            
            if mod(order,2)==1
                maxIndex = (order+1)/2;
            else
                maxIndex = order/2;
            end
            
            O = cell(maxIndex+1,maxIndex+1);
            
            O_11_unnormalised = expm( obj.H_1_tensor );
            normalisation_factor = trace(O_11_unnormalised);
            
            O{1,1} = reshape(  O_11_unnormalised/normalisation_factor , [1,d,d,1] ) ;
            
            
            
            %N=number of free bonds
            for N=1:order
                if mod(N,2)==1
                    current_max_index = (N-1)/2;  %used to contract the tensor
                    double_update(N,current_max_index);
                    
                else
                    current_max_index = N/2;
                    single_update(N,current_max_index)
                end
                
            end
            
            
            
            MPO = mpo_cell_2_matrix(O,maxIndex,d);
            
            
            function  double_update(N,current_max_index)
                %step similar to
                % 0--|--1--|--2--|--1--|--0 = exp(H_12+H_23+H_34) - (0--|--|--|--|--0)
                
                RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1))- contract_O(N,O,current_max_index,d);
                RHS_Matrix = reshape( permute(RHS_Tensor, site_ordering_permute(N+1) ),...
                    dimension_vector(d^2,N+1));%group per ij index
                
                %search x st
                % left*x = res
                % y*right = x
                
                M=current_max_index;
                
                if M ~=0
                    
                    left_list = cell( 1,M );
                    right_list = cell( 1,M );
                    contract_list = cell( 1,M );
                    
                    for i = 1:M
                        left_list{i} = O{i,i+1};
                        right_list{end-i+1} = O{i+1,i};
                        
                        contract_list{i} = [i,-(2*i),-(2*i+1),i+1];
                    end
                    
                    contract_list{1}(1)=-1;
                    contract_list{end}(4)=-(2*M+2);
                    
                    
                    left = reshape( ncon(left_list,contract_list) , [d^(2*M),d^(2*M)]);
                    right = reshape( ncon(right_list,contract_list) , [d^(2*M),d^(2*M)]);
                    res = reshape( RHS_Matrix, [d^(2*M),d^(2* (N+1-M))]);
                    
                    x = left\res;
                    x = reshape(x, [d^(2* (N+1-M)),d^(2*M)] );
                    y = x/right;
                    
                    
                    
                    
                    
                    new_parts = reshape(y, [d^(2*M),d^2,d^2,d^(2*M)]);
                else
                    new_parts = reshape(RHS_Matrix, [d^(2*M),d^2,d^2,d^(2*M)] );
                end
                
                new_parts_sym = reshape( permute( new_parts , [1,2,4,3]), [d^(2*M+2),d^(2*M+2)]);
                new_parts_sym = 0.5*( new_parts_sym+new_parts_sym');
                
                
                
                switch p.Results.method
                    case "diag"
                        %hermitian so eigs exists
                        [U, S] = eig(new_parts_sym, 'vector');
                        
                        [S_2,ind] = sort(S,'descend','ComparisonMethod','abs');
                        U2= U(:,ind);
                        
                        fprintf("S11 %.4e Snn %.4e ",abs(S_2(1)),abs(S_2(2*M+2)));
                        
                        P = sign(S_2);
                        
                        %tol = 1e-14;
                        %S_2(S_2<0)=0;
                        
                        %S_2
                        
                        S_sqrt = sqrt( S_2  );
                        
                        
                        left = U2*diag(S_sqrt);
                        right = diag(P.*S_sqrt)*U2';
                    case "svd"
                        [U,S,V] = svd(new_parts_sym);
                        %S
                        sqrt_S = S.^(0.5);
                        
                        left = U*sqrt_S;
                        right = sqrt_S*V';
                        
                        P = eye(2*M+2);
                    otherwise
                        
                        
                end
                
                O_l =reshape(left, [d^(2*M),d,d,d^(2*(M+1))]);
                O_r =permute( reshape( right, [d^(2*M+2),d^(2*M),d,d])  , [1,3,4,2]);
                
                %
                %                  [U,S,V] = svd(new_parts);
                %                 %to low eigenvalues lead to bad inverses. change
                %                 %normalisation factor if appropriate.
                %
                %                 fprintf("N=%d S(1,1) %.4e S(n,n)  %.4e",N, S(1,1), S(d^(2*(M+1)),d^(2*(M+1))) );
                %
                %
                %                 min = 1e3;
                %                 max = 1e4;
                %                    if S(1,1) < min || S(1,1)>max
                %                         %fprintf("Sii %.2e",S(i,i)) ;
                %                         [O,normalisation_factor] = change_normalisation( O,normalisation_factor, ( (max-min)/(2*S(1,1)) )^(1/(N+1))   );
                %                         double_update(N,current_max_index);
                %                         return
                %
                %                    end
                %               %end
                %
                % %                for i = 1: d^(2*(M+1))
                % %                    if S(i,i) < 1e-7
                % %                         %fprintf("Sii %.2e",S(i,i)) ;
                % %                         %[O,normalisation_factor] = change_normalisation( O,normalisation_factor, (1/S(1,1))^(1/(N+1))   );
                % %
                % %                         S(i,i) = 1e-7;
                % %
                % %                    end
                % %               end
                %
                %               sqrt_X = U*(S.^(1/2))*V';
                %
                %
                %                 O_l =reshape( sqrt_X, [d^(2*M),d,d,d^(2*(M+1))]);
                %                 O_r =  reshape( sqrt_X', [d^(2*(M+1)),d,d,d^(2*M)]);
                
                %check if error is not larger than the good.
                
                %                 if M ~= 0
                %                     %test_y = reshape(y,[d^(2*M),d^4,d^(2*M)]);
                %                     y_err =  ncon(  {left,O_l,O_r,right}, {[-1,1], [1,-2,-3,2],[2,-4,-5,3],[3,-6]})-reshape( RHS_Matrix ,[d^(2*M),d,d,d,d,d^(2*M)]);
                %                     fprintf("sym err %e",tensor_norm(y_err));
                %
                %
                %
                %                 end
                
                
                O{M+1,M+2} = O_l;
                O{M+2,M+1} =O_r;
                
            end
            
            
            function single_update(N,current_max_index)
%                 %             %step 4:
%                 % 0--|--1--|--2--|--2--|--1--|--0 = exp(H_12+H_23+H_34+H_45) - (0--|--|--|--|--|--0)
%                 
%                 
%                 
%                 RHS_Tensor = H_exp(obj,N+1)/(normalisation_factor^(N+2)) - contract_O(N+1,O,current_max_index,d);
%                 RHS_Matrix = reshape( permute(RHS_Tensor, site_ordering_permute(N+2)),...
%                     dimension_vector(d^2,N+2)); %group per ij index
%                 
%                 %search x st
%                 % left*x = res
%                 % y*right = x
%                 M=current_max_index;
%                 
%                 left_list = cell( 1,M );
%                 right_list = cell( 1,M );
%                 contract_list = cell( 1,M );
%                 
%                 for i = 1:M
%                     left_list{i} = O{i,i+1};
%                     right_list{end-i+1} = O{i+1,i};
%                     
%                     contract_list{i} = [i,-(2*i),-(2*i+1),i+1];
%                 end
%                 
%                 contract_list{1}(1)=-1;
%                 contract_list{end}(4)=-(2*M+2);
%                 
%                 
%                 left = reshape( ncon(left_list,contract_list) , [d^(2*M),d^(2*M)]);
%                 right = reshape( ncon(right_list,contract_list) , [d^(2*M),d^(2*M)]);
%                 res = reshape( RHS_Matrix, [d^(2*M),d^(2*M+4  )]);
%                 
%                 x = left\res;
%                 x = reshape(x, [d^(2*M+4),d^(2*M)]);
%                 y = x/right;
%                 
%                 new_part = reshape(y,  [d^(2*M+2 ),d^(2*M+2)]);
%                 %new_part_perm = permute(new_part, [1,2,4,3]);
%                 %temp = reshape( new_part_perm, [d^(2*M+2),d^(2*M+2)]);
%                 
%                 %herm = 0.5*(temp+temp');
%                 
%                 [U, S] = eig(new_part, 'vector');
%                 
%                 [S_2,ind] = sort(S,'descend');
%                 U2= U(:,ind);
%                 
%                 %fprintf("S11 %.4e Snn %.4e ",abs(S_2(1)),abs(S_2(2*M+2)));
%                 
%                 %P = sign(S_2);
%                 
%                 %tol = 1e-14;
%                 %S_2(S_2<0)=0;
%                 
%                 %S_2
%                 
%                 S_sqrt = sqrt( S_2  );
%                 
%                 
%                 U_3 = U2(:,1:2*M);
%                 
%                 
%                 %fprintf( "orig_err %.4e new_err %.4e ",orig_err,new_err );
%                 
%                 O{M+1,M+1} = reshape( U2*S_sqrt*U3',[d^(2*M),d,d,d^(2*M)])   ;
%                 
%                 %                 if new_err > orig_err
%                 %                     fprintf("bad mpo order %d, ",N);
%                 %                     %return
%                 %                 else
%                 %
%                 %                 end
%                 
%                 
                
                RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1)) - contract_O(N,O,current_max_index,d);
                RHS_Matrix = reshape( permute(RHS_Tensor, site_ordering_permute(N+1)),...
                                        dimension_vector(d^2,N+1)); %group per ij index

                %search x st 
                % left*x = res
                % y*right = x
                 M=current_max_index;
                
                left_list = cell( 1,M );
                right_list = cell( 1,M );
                contract_list = cell( 1,M );
                
                for i = 1:M
                    left_list{i} = O{i,i+1};
                    right_list{end-i+1} = O{i+1,i};
                    
                    contract_list{i} = [i,-(2*i),-(2*i+1),i+1];
                end
                
                contract_list{1}(1)=-1;
                contract_list{end}(4)=-(2*M+2);
                

                left = reshape( ncon(left_list,contract_list) , [d^(2*M),d^(2*M)]);
                right = reshape( ncon(right_list,contract_list) , [d^(2*M),d^(2*M)]);
                res = reshape( RHS_Matrix, [d^(2*M),d^(2*M+2)]);

                x = left\res;
                x = reshape(x, [d^(2*M+2),d^(2*M)]);
                y = x/right;                                   

% 
%                 err = ncon( { left,reshape(y, [d^(2*M),d^2,d^(2*M)]),right},{[-1,1],[1,-2,2],[2,-3]})- reshape(res, [d^(2*M),d^2,d^(2*M)]);
%                 
%                 fprintf("assym err %e",tensor_norm(err));
%                 
                O{M+1,M+1} = reshape(y, [d^(2*M),d,d,d^(2*M)]);  
                
                
            end
            
            
            function [O,normalisation_factor] = change_normalisation(  O,normalisation_factor,change)
                max_dim = size(O,1);
                
                for i =1:max_dim
                    
                    for j = 1:max_dim
                        O{i,j} = O{i,j}*change;
                    end
                end
                
                normalisation_factor = normalisation_factor/change;
                
            end
        end
        
        function [normalisation_factor,MPO] = type_02(obj, testing)
            % Make a cell from O. It holds the tensor elements
            % every entry is 4d nxdxdxm with n and m the bond dimension for the
            % corresponding bond. dimension x1 at the end not shown by matlab
            %this type generates no 1--|--1 and 2--|--2 blocks
            d = obj.dim;
            
            maxIndex = 4;
            O = cell(maxIndex+1,maxIndex+1);
            
            O_11_unnormalised = expm( obj.H_1_tensor );
            normalisation_factor = trace(O_11_unnormalised);
            
            O{1,1} = reshape(   O_11_unnormalised/normalisation_factor, [1,d,d,1] ) ;
            
            %step 1:
            % 0--|--1--|--0 = exp(H_12) - (0--|--|--0 )
            N = 1;                  %number of free bonds
            current_max_index = 0;  %used to contract the tensor
            
            RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1))- contract_O(N, O ,current_max_index,d);
            RHS_Matrix_site = reshape( permute( RHS_Tensor , site_ordering_permute(N+1) ),...
                [d^2,d^2] ); %ready to svd
            
            [U,S,V] = svd(RHS_Matrix_site);
            
            a_S = average(S);
            
            O{0+1,1+1} = reshape( U* a_S , [1,d,d,d^2]);
            O{1+1,0+1} = reshape( S/a_S * V', [d^2,d,d,1]);
            
            if testing==1
                err = tensor_norm( ncon( {O{0+1,1+1},O{1+1,0+1}}, {[-1,-2,-4,1],[1,-3,-5,-6]}, [1])-RHS_Tensor);
                fprintf("err 01 = %d\n",err);
            end
            
            %step 2 :
            % 0--|--1--|--2--|--0 = exp(H12+H23)- (0--|--|--|--0)
            N = 2;                  %number of free bonds
            current_max_index = 1;  %used to contract the tensor
            
            RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1)) - contract_O(N,O,current_max_index,d);
            RHS_Matrix_site = reshape(permute(RHS_Tensor, site_ordering_permute(N+1) ),...
                [d^2,d^2,d^2] );
            
            left = reshape( O{1,2}, [d^2,d^2]);
            right = reshape(RHS_Matrix_site, [d^2,d^4]);
            
            %https://nl.mathworks.com/help/matlab/ref/mldivide.html
            x = left\right;
            
            new_part = reshape( x, [d^4,d^2] );
            
            [U,S,V] = svd(new_part);
            a_S = average(S);
            
            O{1+1,2+1} = reshape( U*a_S, [d^2,d,d,d^4]);
            O{2+1,0+1} = reshape( S/a_S*V', [d^4,d,d,1]);
            
            if testing==1
                err = tensor_norm( ncon( { O{1,2}, O{2,3}, O{3,1} }, {[-1,-2,-5,1],[1,-3,-6,2],[2,-4,-7,-8]})...
                    -RHS_Tensor);
                fprintf("err 0120 = %d\n",err);
            end
            
            %step 3 :
            % 0--|--1--|--2--|--3--|--0 = exp(H12+H23+h34)- (0--|--|--|--|--0)
            N = 3;                  %number of free bonds
            current_max_index = 2;  %used to contract the tensor
            
            RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1)) - contract_O(N,O,current_max_index,d);
            RHS_Matrix_site = reshape(permute(RHS_Tensor, site_ordering_permute(N+1) ),...
                [d^2,d^2,d^2,d^2] );
            
            left = reshape( ncon( {O{1,2},O{2,3}},{[-1,-2,-3,1],[1,-4,-5,-6]} ),...
                [d^4,d^4]);
            right = reshape(RHS_Matrix_site, [d^4,d^4]);
            
            %https://nl.mathworks.com/help/matlab/ref/mldivide.html
            x = left\right;
            
            new_part =  reshape( x, [d^6,d^2] );
            
            [U,S,V] = svd( new_part);
            a_S = average(S);
            
            O{2+1,3+1} = reshape( U*a_S, [d^4,d,d,d^6]);
            O{3+1,0+1} = reshape( S/a_S*V', [d^6,d,d,1]);
            
            if testing==1
                err = tensor_norm( ncon( { O{1,2}, O{2,3}, O{3,4},O{4,1} }, {[-1,-2,-6,1],[1,-3,-7,2],[2,-4,-8,3] ,[3,-5,-9 ,-10 ]} )...
                    -RHS_Tensor);
                fprintf("err 01230 = %d\n",err);
            end
            
            %step 4 :
            % 0--|--1--|--2--|--3--|--4--|--0 = exp(H12+H23+h34)- (0--|--|--|--|--|--0)
            N = 4;                  %number of free bonds
            current_max_index = 3;  %used to contract the tensor
            
            RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1)) - ...
                contract_O(N,O,current_max_index,d);
            RHS_Matrix_site = reshape(permute(RHS_Tensor, site_ordering_permute(N+1) ),...
                [d^2,d^2,d^2,d^2,d^2] );
            
            left = reshape( ncon( {O{1,2},O{2,3},O{3,4}},{[-1,-2,-3,1],[1,-4,-5,2],[2,-6,-7,-8]} ),...
                [d^6,d^6]);
            right = reshape(RHS_Matrix_site, [d^6,d^4]);
            
            %https://nl.mathworks.com/help/matlab/ref/mldivide.html
            x = left\right;
            
            new_part =  reshape( x , [d^8,d^2] );
            
            
            [U,S,V] = svd( new_part);
            a_S = average(S);
            
            O{3+1,4+1} = reshape( U*a_S, [d^6,d,d,d^8]);
            O{4+1,0+1} = reshape( S/a_S*V', [d^8,d,d,1]);
            
            if testing==1
                err = tensor_norm( ncon( { O{1,2}, O{2,3}, O{3,4},O{4,5},O{5,1} }, {[-1,-2,-7,1],[1,-3,-8,2],[2,-4,-9,3] ,[3,-5,-10 ,4 ],[4,-6,-11,-12]} )...
                    -RHS_Tensor);
                fprintf("err 012340 = %d\n",err);
            end
            
            
            MPO = mpo_cell_2_matrix(O,maxIndex,d);
            
        end
        
        
        function [normalisation_factor,MPO] = type_03(obj,order ,testing)
            % Make a cell from O. It holds the tensor elements
            % every entry is 4d nxdxdxm with n and m the bond dimension for the
            % corresponding bond. dimension x1 at the end not shown by matlab
            %this type generates no 1--|--1 and 2--|--2 blocks
            %uses virtual levels:
            %0--|--1--|--0
            %0--|--1'--|--2'--|--0
            %0--|--1''--|--2''--|--3''--|--0
            
            
            %representation blocks in block matrix with dims
            
            %
            %       1         d^2     d^2    d^2      d^2     d^4         d^2
            %
            %  1     00        01'    |01''   0       |01'''   0           0       |
            %  d^2   1'0       0      |0      0       |0       0           0       |
            %       __________________|               |                            |
            %  d^2   0         0       0       1''2'' |0       0           0       |
            %  d^2   2''0      0       0       0      |0       0           0       |
            %       __________________________________|0       0           0       |
            %  d^2   0         0       0       0       0       1'''2'''    0       |
            %  d^4   0         0       0       0       0       0           2'''3'''|
            %  d^2   3'''0     0       0       0       0       0           0       |
            %       _______________________________________________________________|
            
            
            
            
            d = obj.dim;
            
            total_dim=1;
            
            for k = 1:order
                for i =1:k
                    total_dim =  total_dim +internal_dim(i,k,d);
                end
            end
            
            
            %total_dim = 1+ d^2 + (d^2+d^2)+ (d^2+d^4+d^2)+(d^2+d^4+d^4+d^2);
            
            
            left_vect= zeros(1,total_dim);
            left_vect(1)=1;
            right_vect= zeros(total_dim,1);
            right_vect(1)=1;
            
            MPO = zeros(total_dim,d,d,total_dim);
            
            % 0
            unnorm = expm(obj.H_1_tensor);
            normalisation_factor = trace(unnorm);
            
            O_00_0 = reshape(unnorm/normalisation_factor, [1,d,d,1] );
            
            MPO = add_block_to_tensor(MPO,0,0,0, O_00_0 ,d  );
            
            
            %step 1:
            % 0--|--1'--|--0 = exp(H_12) - (0--|--|--0 )
            N = 1;                  %number accents = number of free bonds
            
            
            RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1))- contract_MPO(MPO,N,left_vect,right_vect);
            RHS_Matrix_site = reshape( permute( RHS_Tensor , site_ordering_permute(N+1) ),...
                [d^2,d^2] ); %ready to svd
            
            O_01_1=  reshape(eye(d^2),[1,d,d,d^2]);
            O_10_1=  reshape(RHS_Matrix_site, [d^2,d,d,1]);
            
            MPO = add_block_to_tensor(MPO,0,1,N, O_01_1, d  );
            MPO = add_block_to_tensor(MPO,1,0,N, O_10_1, d  );
            
            
            if testing==1
                err = tensor_norm( ncon( {O_01_1,O_10_1}, {[-1,-2,-4,1],[1,-3,-5,-6]})-RHS_Tensor);
                fprintf("err 01 = %d\n",err);
            end
            
            %all other orders
            for N=2:order
                construct_level(N)
            end
            
            
            function construct_level(N)
                %step 4 :
                % 0--|--1''''--|--2''''--|--3''''--|--4''''--|--0 = exp(H12+H23+h34)- (0--|--|--|--|--|--0)
                %N=4
                
                %N=3
                %0--|--1'''--|--2'''--|--3'''---|---0
                
                if mod(N,2)==0
                    M_l = N/2;
                    M_r= N/2;
                else
                    M_l = (N-1)/2;
                    M_r= (N-1)/2+1;
                end
                
                RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1)) - contract_MPO(MPO,N,left_vect,right_vect);
                RHS_Matrix_site = reshape(permute(RHS_Tensor, site_ordering_permute(N+1) ),...
                    [d^(2*M_l),d,d,d^(2*M_r) ] );
                
                for k = 0:M_l-1
                    MPO = add_block_to_tensor(MPO,k,k+1,N, reshape( eye(d^(2*(k+1)) ), [d^(2*k),d,d,d^(2*(k+1))]) , d  );
                end
                
                MPO = add_block_to_tensor(MPO,M_l,M_l+1,N, reshape( RHS_Matrix_site, [d^(2*M_l) ,d,d,d^(2*M_r)]) , d  );
                
                for k = M_l+1:M_l+M_r-1
                    MPO = add_block_to_tensor(MPO,k,k+1,N, reshape( eye(d^(2*(N-k+1) ) ), [d^(2*( N-k+1)),d,d,d^(2*(N-k))]) , d  );
                end
                
                MPO = add_block_to_tensor(MPO,N,0,N, reshape( eye(d^2) , [d^2,d,d,1] ) , d  );
            end
            
            function y=  internal_dim(i,k,d)
                y=d^( 2* min(i,k-i+1)  );
            end
            
            function T = contract_MPO(MPO,N,left,right)
                M=N+3;
                
                tensors = cell(1,M);
                tensors(2:end-1)= {MPO};
                tensors{1}=left;
                tensors{M} = right;
                
                
                leg_list = cell(1,M);
                leg_list{1} = [-1,1];
                leg_list{M} = [M-1,-2*(M-1)];
                
                for i = 2:M-1
                    leg_list{i} = [i-1,-(i), -(i+M-2)  ,i  ] ;
                end
                
                
                T = ncon( tensors,leg_list   );
                
            end
            
            
            
            %test function for type_03
            function test_add_block_to_tensor
                d=2;
                dim = 1+d^2+(d^2+d^2) + (d^2+d^4+d^2);
                O = zeros(dim,d,d,dim);
                
                b_00_1 = ones(1,d,d,1)*0.5;
                O = add_block_to_tensor(O,0,0,0,b_00_1,d);
                
                b_01_1 = ones(1,d,d,d^2);
                O = add_block_to_tensor(O,0,1,1,b_01_1,d);
                
                b_10_1 = ones(d^2,d,d,1)*2;
                O = add_block_to_tensor(O,1,0,1,b_10_1,d);
                
                b_01_2 = ones(1,d,d,d^2)*3;
                O = add_block_to_tensor(O,0,1,2,b_01_2,d);
                
                b_12_2 = ones(d^2,d,d,d^2)*5;
                O = add_block_to_tensor(O,1,2,2,b_12_2,d);
                
                b_20_2 = ones(d^2,d,d,1)*4;
                O = add_block_to_tensor(O,2,0,2,b_20_2,d);
                
                
                b_01_3 = ones(1,d,d,d^2)*6;
                O = add_block_to_tensor(O,0,1,3,b_01_3,d);
                
                b_12_3 = ones(d^2,d,d,d^4)*7;
                O = add_block_to_tensor(O,1,2,3,b_12_3,d);
                
                b_23_3 = ones(d^4,d,d,d^2)*8;
                O = add_block_to_tensor(O,2,3,3,b_23_3,d);
                
                b_30_3 = ones(d^2,d,d,1)*9;
                O = add_block_to_tensor(O,3,0,3,b_30_3,d);
                
                Z= reshape(O(:,1,1,:),[dim,dim]);
                
            end
            
            
            %representation blocks in block matrix.
            
            %
            %        1         d^2     d^2    d^2      d^2     d^4         d^2
            %
            %  1     00        01'    |01''   0       |01'''   0           0       |
            %  d^2   1'0       0      |0      0       |0       0           0       |
            %       __________________|               |                            |
            %  d^2   0         0       0       1''2'' |0       0           0       |
            %  d^2   2''0      0       0       0      |0       0           0       |
            %       __________________________________|0       0           0       |
            %  d^2   0         0       0       0       0       1'''2'''    0       |
            %  d^4   0         0       0       0       0       0           2'''3'''|
            %  d^2   3'''0     0       0       0       0       0           0       |
            %       _______________________________________________________________|
            
            
            %add block  i(k)--|--j(k) to the matrix
            %O (dim,d,d,dim)
            function O = add_block_to_tensor(O,i,j,k,block,d  )
                
                block_start =1;%00 block
                
                if i==0 && j==0
                    O(1,:,:,1) = block(:,:,:,:);
                    return;
                end
                
                %todo make explicit formula
                for  s= 1:k-1
                    for t=1:s
                        block_start= block_start+ internal_dim(t,s,d );
                    end
                end
                
                if i==0
                    y= block_start+1;
                    x=1;
                    
                else
                    if j==0
                        x = block_start+1;
                        for s = 1:i-1
                            x = x+internal_dim(s,k,d);
                        end
                        y= 1;
                    else
                        x = block_start+1;
                        for s = 1:i-1
                            x = x+internal_dim(s,k,d);
                        end
                        
                        y = block_start+1;
                        for t = 1:j-1
                            y = y+internal_dim(t,k,d);
                        end
                    end
                end
                
                xdim = size(block,1);
                ydim = size(block,4);
                O(x:x+xdim-1,:,:,y:y+ydim-1) = O(x:x+xdim-1,:,:, y:y+ydim-1)+block(:,:,:,:);
                
                % k is subspace dim
                function y=  internal_dim(i,k,d)
                    y=d^( 2* min(i,k-i+1)  );
                end
            end
            
            
            
        end
        
        
        function [normalisation_factor,MPO] = type_05(obj,order,opt)
            %this mpo limits the possibilities and works with unitary transformations
            %
            %  unitary: (S_(n n+1)^(i alpha j)_beta) ' = T_(n+1 n)^ alpha _ (i beta j)
            % matrix is constructed st only these combinations can happen
            %  S01--S12--...--S(n-1 n)--D_nn--n--S(n-1 n)'--...--S21'--S10'    with D a real diagonal matrix
            %  S01--S12--...--S(n-1 n)--(T_nn--)^m--S(n-1 n)'--...--S21'--S10' with
            %  T_nn a tesnsor
            
            % cosntruction matrix,  ' means hermitian conjugate. Example for blocks uo
            % until 2, can be generalised for any n
            %      B 0 |   B1     |     B2     |     B3    |    B4
            
            % dims
            %        1   d^2  d^4   d^2   d^4    d^2  d^4    d^2      d^4
            %
            % 1    T_00| S01  0   |-2*S01  0   | S01  0    | S01*D11   0
            %      ____|          |            | 0    S12  | 0         0
            % d^2  S01'  0    S12 | 0 -2*S02   | 0    0    | 0         0
            % d^4  0     S12' 0   | 0      0   | 0    0    | 0         0
            %      _______________| 0      0   | 0    0    | 0         0
            % d^2  S01'  0    0     0      0   | 0    0    | 0         0
            % d^4  0     S02' 0     0      0   | 0    0    | 0         0
            %      ____________________________| 0    0    | 0         0
            % d^2  S01'  0    0     0      0     T11  0    | 0         0
            % d^4  0     S12' 0     0      0     0    T22  | 0         0
            %      ________________________________________| 0         0
            % d^2  S01'  0    0     0      0     0    0      0         1/D11*S12*D22
            % d^4  0     0    0     0      0     0    0      S12'      0
            
            %intuition matrix: B4: create blocks with D,  B3 with B1: create T_nn
            %sequences. B1 is is used to create S(n n+1) sequence beforen T_nn^m
            % block B2 corrects for spurious multiplications in B1 and B3 that do not
            % involve any T_nn
            
            %setup
            d=obj.dim;
               
            if mod(order,2)==1
                max_index = (order+1)/2;
            else
                max_index = order/2;
            end
            
            total_dim =  1+ 4* (  ((d^2)^(max_index+1)-1)/(d^2-1) - 1)  ;
            
            left_vect= zeros(1,total_dim);
            left_vect(1)=1;
            right_vect= zeros(total_dim,1);
            right_vect(1)=1;
            
            MPO = zeros(total_dim,d,d,total_dim);
            
            U_cell= cell(1,max_index) ;
            V_cell= cell(1,max_index) ;
            
            sqrt_Dn = 1;
            
            % 00 block
            N=0;
            
            unnorm = expm(obj.H_1_tensor);
            normalisation_factor = trace(unnorm);
            
            T00 = reshape(unnorm/normalisation_factor, [1,d,d,1] );
            MPO = add_block_Tn(MPO,N,T00,order,d );
            
            %other blocks
            %N=number of free bonds
            for N=1:order
                if mod(N,2)==1
                    current_max_index = (N-1)/2;  %used to contract the tensor
                    [sqrt_Dn]=double_update(current_max_index,N,sqrt_Dn,max_index);           
                else
                    current_max_index = N/2;
                    single_update(current_max_index,N,max_index)
                end
                
            end
            
            Z=reshape(  MPO(:,1,1,:), [total_dim,total_dim ]);
            
            %MPO = mpo_cell_2_matrix(O,maxIndex,d);
            
            %n current num of bonds
            %N max free bonds (order)
            function [sqrt_Dn] = double_update(current_max_index,N,sqrt_Dnm,max_index)
               %  S01--S12--...--S(n-1 n)--D_nn--n--S(n-1 n)'--...--S21'--S10' with D a real diagonal matrix
            
                M=current_max_index;

                RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1))- contract_MPO(MPO,N,left_vect,right_vect);
                RHS_Matrix = reshape( permute(RHS_Tensor, site_ordering_permute(N+1) ),  [d^(2*M),d^2,d^2,d^(2*M) ]  );
                
                res = RHS_Matrix;
                
                for i = 1:current_max_index
                     %apply one at a time
                    U= U_cell{i};
                    V= V_cell{i};
                    res = ncon(  {U',  reshape(res, [ d^(2*i )  , d^(2*M-2*i+2) , d^(2*M-2*i+2) , d^(2*i ) ]) , V}, { [-1,1],[1,-2,-3,2],[2,-4]});  
    
                end

                
                new_parts = reshape( res ,[d^(2*M+2),d^(2*M+2)]);
                %new_parts_sym = 0.5*( new_parts_sym+new_parts_sym');%todo check wheter this is good
                
                [U,Dn,V] = svd(new_parts);
              
                sqrt_Dn = Dn.^(1/2);
                
                
              
                    
                Sn = reshape(U,[d^(2*M),d,d,d^(2*M+2)]);
                SnD =  reshape(V',[d^(2*M+2),d,d,d^(2*M)]);
                

               
                U_cell{1+current_max_index}= reshape( U, [d^(2*M+2),d^(2*M+2)]);
                V_cell{1+current_max_index}= reshape( V, [d^(2*M+2),d^(2*M+2)]);
                
                
                %undo previous Sn
   

                MPO = add_block_05(MPO,current_max_index+1,Sn,SnD,sqrt_Dn,sqrt_Dnm,max_index,d );
               
                %for debugging purposes
                
                
            end
            
            
            function single_update(current_max_index,N,max_index)
                %  S01--S12--...--S(n-1 n)--(T_nn--)^m--S(n n-1)'--...--S21'--S10' with
               
                RHS_Tensor = H_exp(obj,N)/(normalisation_factor^(N+1))...
                    - contract_MPO(MPO,N,left_vect,right_vect);
                
                
                RHS_Matrix = reshape( permute(RHS_Tensor, site_ordering_permute(N+1)),...
                    dimension_vector(d^2,N+1)); %group per ij index
                
                M=current_max_index;
                
              
                res= RHS_Matrix;
                
                for i = 1:current_max_index
                     %apply one at a time
                    U= U_cell{i};
                    V= V_cell{i};
                    res = ncon(  {U',  reshape(res, [ d^(2*i )  ,d^(2*M+1-2*i),d^(2*M+1-2*i), d^(2*i ) ]) , V}, { [-1,1],[1,-2,-3,2],[2,-4]});  
    
                end
                

                              
                %add_block_Tn(MPO,n,Tn,max_index,d )
                MPO = add_block_Tn(MPO,current_max_index,res,max_index,d );
                
                %for debugging purposes
                Z=reshape(  MPO(:,1,1,:), [total_dim,total_dim ]);
                
            end                      
           
            
            
           
            
            function T = contract_MPO(MPO,N,left,right)
                M=N+3;
                
                tensors = cell(1,M);
                tensors(2:end-1)= {MPO};
                tensors{1}=left;
                tensors{M} = right;
                
                
                leg_list = cell(1,M);
                leg_list{1} = [-1,1];
                leg_list{M} = [M-1,-2*(M-1)];
                
                for i = 2:M-1
                    leg_list{i} = [i-1,-(i), -(i+M-2)  ,i  ] ;
                end
                
                T = ncon( tensors,leg_list   );
            end
             
           
            %add blocks to mpo like in discription
            function MPO =  add_block_05(MPO,n,Sn,SnD,sqrt_Dn,sqrt_Dnm,max_index,d )
                
%                 if n==0
%                     MPO(1,:,:,1)=Tn;
%                     return
%                 end
                
                %horizontal
                block_start_x = 1;
                
                %B=1 case
                B=1;
                block_start_y = get_B_start(B,d,max_index);
                internal_diag = geom_sum(n-1,d);
                internal_diag_m = geom_sum(n-2,d)+1;
                MPO(1+ block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2*n)-1)  ,:,:, block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2*n-2)-1) ) = SnD;
                %mirror case for
                MPO(block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2*n-2)-1),:,:,1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2*n)-1) ) = Sn;
                
%                 if n==1
                    %B=2 case
                    B=2;
                    block_start_y = get_B_start(B,d,max_index);
                    internal_diag = geom_sum(n-1,d);
                    internal_diag_m = geom_sum(n-2,d)+1;
                    MPO( 1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2*n)-1),:,:, block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2*n-2)-1) ) = SnD;
                    %mirror case for
                    MPO(block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2*n-2)-1),:,:,1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2*n)-1) ) = -2* Sn;
%                 else
%                     B=2;
%                     block_start_y=get_B_start(B,d,max_index);
%                     internal_diag = geom_sum(n-1,d);
%                     internal_diag_m = geom_sum(n-2,d);
%                     
%                     
%                     MPO(1+block_start_y+internal_diag_m :1+block_start_y+internal_diag_m + (d^(2*n-2)-1),:,:,1+block_start_y+internal_diag:1+block_start_y+internal_diag+ (d^(2*n)-1) ) = Sn;
%                     MPO(1+block_start_y+internal_diag:1+block_start_y+internal_diag+ (d^(2*n)-1),:,:, 1+block_start_y+internal_diag_m :1+block_start_y+internal_diag_m + (d^(2*n-2)-1) ) = SnD;
%   
%                 end
                %B=3 case
                B=3;
                block_start_y = get_B_start(B,d,max_index);
                internal_diag = geom_sum(n-1,d);
                internal_diag_m = geom_sum(n-2,d)+1;
                MPO(1+ block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2*n)-1),:,:,block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2*n-2)-1) ) = SnD;
                %mirror case for
                MPO( block_start_x+internal_diag_m:block_start_x+internal_diag_m+(d^(2*n-2)-1),:,:,1+block_start_y+internal_diag:1+block_start_y+internal_diag+(d^(2*n)-1) ) = Sn;
                
%                 z = block_start_y+internal_diag;
%                 MPO(1+z:1+z+(d^(2*n)-1) ,:,:, 1+z:1+z+(d^(2*n)-1) ) = Tn;
%                 
                
                if n ==1
                    
                    SD = ncon( {Sn, sqrt_Dn},{ [-1,-2,-3,1],[1,-4]});
                    
                    SDn = ncon( {sqrt_Dn,SnD },{ [-1,1],[1,-2,-3,-4]});
                    
                    block_start_y=get_B_start(4,d,max_index);
                    %MPO
                    MPO(1,:,:,1+block_start_y:1+block_start_y+ (d^2-1) ) = SD;
                    MPO(1+block_start_y:1+block_start_y+ (d^2-1),:,:,1 ) = SDn;
                    
                else  %todo
                    block_start_y=get_B_start(4,d,max_index);
                    internal_diag = geom_sum(n-1,d);
                    internal_diag_m = geom_sum(n-2,d)+1;
                    
                    
                    sqrt_Dnm_inv = sqrt_Dnm^-1; %is diagonal
                    
                    SD = ncon( {sqrt_Dnm_inv, Sn, sqrt_Dn},{ [-1,1] [1,-2,-3,2],[2,-4]});
                    
                    SDn =  ncon( {sqrt_Dn, SnD, sqrt_Dnm_inv},{ [-1,1] [1,-2,-3,2],[2,-4]});
                    
                    MPO(block_start_y+internal_diag_m:block_start_y+internal_diag_m+ (d^(2*n-2)-1),:,:,1+block_start_y+internal_diag:1+block_start_y+internal_diag+ (d^(2*n)-1) ) = SD;
                    MPO(1+block_start_y+internal_diag:1+block_start_y+internal_diag+ (d^(2*n)-1),:,:, block_start_y+internal_diag_m:block_start_y+internal_diag_m+ (d^(2*n-2)-1) ) = SDn;
                end
                
                function y= get_B_start(B,d,N)
                    block_dim = geom_sum(N,d);

                    y= 1 + (B-1)*block_dim;
                end

                function y= geom_sum(N,d)
                    y=  ((d^2)^(N+1)-1)/(d^2-1) - 1;
                end
            end
                     
            
            
            function MPO =  add_block_Tn(MPO,n,Tn,max_index,d )
                
                if n==0
                    MPO(1,:,:,1)=Tn;
                    return
                end
              
                B=3;
                block_start_y = get_B_start(B,d,max_index);
                internal_diag = geom_sum(n-1,d);
               
                z = block_start_y+internal_diag;
                MPO(1+z:1+z+(d^(2*n)-1) ,:,:, 1+z:1+z+(d^(2*n)-1) ) = Tn;
                
                function y= get_B_start(B,d,N)
                    block_dim = geom_sum(N,d);
                    
                    y= 1 + (B-1)*block_dim;
                end
                
                function y= geom_sum(N,d)
                    y=  ((d^2)^(N+1)-1)/(d^2-1) - 1;
                end
            end    
                
                
            
            %test code for previous function, everything works fine
            %old version where t also was added in same fn
            function test_add_block
                d=2;
                N=3;
                total_dim =  1+ 4* (  ((d^2)^(N+1)-1)/(d^2-1) - 1)  ;
                MPO = zeros( total_dim,d,d,total_dim );
                
                %     add_block_05(MPO,n,Sn,SnD,Tn,            Dn,Dnm,N,d )
                MPO = add_block_05(MPO,0,0, 0, 0.3* ones(1,d,d,1) ,0,0,N,2);
                
                %S1
                N=1;
                S01 = ones(1,d,d,d^2)*0.5454;
                S10 = ones(d^2,d,d,1)*0.5454;
                T11 = ones(d^2,d,d,d^2)*3.3333;
                D11=  eye(d^2 );
                %     add_block_05(MPO,n,Sn,SnD,Tn,Dn,Dnm,N,d )
                MPO = add_block_05(MPO,N,S01,S10,T11,D11,1,N,2);
                
                %S2
                N=2;
                S12 = ones(d^2,d,d,d^4)*0.156165;
                S21 = ones(d^4,d,d,d^2)*0.156165;
                T22 = ones(d^4,d,d,d^4)*0.1686412;
                D22=  eye(d^4 );
                %     add_block_05(MPO,n,Sn,SnD,Tn,Dn,Dnm,N,d )
                MPO = add_block_05(MPO,N,S12,S21,T22,D22,D11,N,2);
                
                
                %S3
                N=3;
                S23 = ones(d^4,d,d,d^6)*0.7777;
                S32 = ones(d^6,d,d,d^4)*0.7777;
                T33 = ones(d^6,d,d,d^6)*0.1686412;
                D33=  eye(d^6 );
                %     add_block_05(MPO,n,Sn,SnD,Tn,Dn,Dnm,N,d )
                MPO = add_block_05(MPO,N,S23,S32,T33,D33,D22,N,2);

                Z=reshape( MPO(:,1,1,:), [total_dim,total_dim]) ;
             end

        end
        

        function H_exp = H_exp(obj,N,matrix)
            if nargin <3
                matrix = 0;
            end
            
            
            % return E(H_1_2+..+H_N-1_N) in normal ordering (dimension d^N+1, basis first
            % upper legs, then lower legs
            % so this makes first the tensor T = H     x I x I
            %                                  + I x H     x I ...
            %
            %                                  + S x I x I x I
            %                                  + I x S x I x I ...
            % with S the single site operator and H the 2 site one
            % then reorders, and exponentiates
            d = obj.dim;
            
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
                
                
                H_i = ncon( tensor_list, leg_list);
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
            
            H_exp_matrix = expm( reshape(H, [d^(N+1),d^(N+1)]));
            
            if matrix == 1
                H_exp =  H_exp_matrix;
                return
            end
            
            H_exp = reshape(  H_exp_matrix , dimension_vector(d,2*(N+1),[1,1] ));
            
            
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
        
        
        T_j = ncon (O_tensors,legs);
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




end

function z =geomSum(x,n)
z= (x^(n+1)-1)/(x-1);
end

%s diagonal
function y = average(S)
d= min(size(S,1),size(S,2));
y=0;
for i = 1:d
    y=y+S(i,i) ;
end
y=y/d;
end



