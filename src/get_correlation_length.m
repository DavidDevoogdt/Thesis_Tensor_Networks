
%%
 function [A,xi_inv] = thermal_correlation_helper(MPO_matrix,O,N1,N2,a,nf)

    
    Traced_MPO = ncon( {MPO_matrix},{[-1,1,1,-2]});
    [V,D,W] = eig(Traced_MPO);
    
    [d,ind] = sort(diag(D));
    
    largest_ind = ind(end);
    
    rho_right= V(:,largest_ind);
    rho_left = W(:,largest_ind)';
    lambda = D(largest_ind,largest_ind)
    
    MPO_matrix_norm = MPO_matrix/lambda;
    traced_O_norm = Traced_MPO/lambda;
    
    n = rho_left*rho_right ;
    
    x_array = N1:N2;
    plot_array = zeros(N2-N1+1,1);
    for i = x_array
        MPO_n = traced_O_norm^(i-1);
        c =1/(n)* ncon( {rho_left,MPO_matrix_norm,O,MPO_n,MPO_matrix_norm,O,rho_right}, {  [-1,1], [1, 2,3,4],[3,2],[4,5], [5,6,7,8],[7,6],[8,-2]} );
        plot_array(i-N1+1) = c;
    end
    
    a_array = x_array*a;
    
   
    
    F = fit(  transpose( a_array) , plot_array, 'exp1');

%       semilogy(a_array,plot_array)
%       hold on
%       semilogy(a_array,F(a_array) )
%       hold off

    C =  coeffvalues(F);
    A = C(1);
    xi_inv = - C(2);

     
 end