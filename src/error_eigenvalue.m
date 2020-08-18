function err = error_eigenvalue(mpo_1_list, mpo_2,d)
    %calculates difference of cycle of N mpo's and finds the larges eigenvector
    %no optimisation is done, just brute force
    
    a= generate_cycle(mpo_1_list,d);
    
    if ~iscell(mpo_2)
        b = mpo_2;
    else
        b =generate_cycle(mpo_2,d);
    end
    
   
    
    matrix=a-b;
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
    
    %todo do this in steps
    temp = ncon( mpo_list,leg_list   );
    
    matrix = reshape( temp, [ d^M,d^M ]);
    
end