function err = calculate_error(a,b)

    N = numel(b)^(0.5) ;
    
    

    b = reshape( b , [N,N]);
    a = reshape( a , [N,N]);

    p = 2;

    [~, S1, ~] = svds(a, 30);

    sum_1 = (sum(diag(S1).^p))^(1 / p);

    [~, S2, ~] = svds(b, 30);

    sum_2 = (sum(diag(S2).^p))^(1 / p);

    err = sum_1 / sum_2;
   
end
