function err = error_eigenvalue(mpo_1_list, mpo_2_list,d)
    %calculates difference of cycle of N mpo's and finds the larges eigenvector
    %no optimisation is done, just brute force
    
    matrix  = generate_cycle(mpo_1_list,d) -  ...
        generate_cycle(mpo_2_list,d);
    err = eigs(matrix,1);
    
end

function matrix = generate_cycle(mpo_list,d)
    M = size(mpo_list,2);
   
    leg_list = cell(1,M);
    for i = 1:M
       leg_list{i} = [i,-(i+1), -(i+M+1)  ,i+1  ] ;
    end

    leg_list{1}(1) = -1;
    leg_list{M}(4) = -(2*M+2);
    
    temp = ncon( mpo_list,leg_list   );
    
    matrix = reshape( temp, [ d^M,d^M ]);
    
end