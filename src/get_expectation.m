%correlator
function m = get_expectation(MPO_matrix,O)
    Traced_MPO = ncon( {MPO_matrix},{[-1,1,1,-2]});
    [V,D,W] = eig(Traced_MPO);
    
    
    [d,ind] = sort(diag(D));
    
    largest_ind = ind(end);
    
    rho_right= V(:,largest_ind);
    rho_left = W(:,largest_ind)';
    lambda = D(largest_ind,largest_ind);
    
    n = rho_left*rho_right ;
    
    
    m =  1/(lambda*n ) *ncon( { rho_left, MPO_matrix, O, rho_right}, { [-1,1], [1,2,3,4],[3,2],[4,-2]});
end





 