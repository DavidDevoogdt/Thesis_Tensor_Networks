
%%
 function [A,xi_inv] = thermal_correlation_helper(MPO_matrix,O,N)

    
    Traced_MPO = ncon( {MPO_matrix},{[-1,1,1,-2]});
    [V,D,W] = eig(Traced_MPO);
    
    [d,ind] = sort(diag(D));
    
    largest_ind = ind(end);
    
    rho_right= V(:,largest_ind);
    rho_left = W(:,largest_ind)';
    lambda = D(largest_ind,largest_ind);
    
    n = rho_left*rho_right ;
    
    x_array = 1:N;
    plot_array = zeros(N,1);
    for i = x_array

        MPO_n = Traced_MPO^(i-1);
        c =1/(n*lambda^(i+1))* ncon( {rho_left,MPO_matrix,O,MPO_n,MPO_matrix,O,rho_right}, {  [-1,1], [1, 2,3,4],[3,2],[4,5], [5,6,7,8],[7,6],[8,-2]} );

        plot_array(i) = c;
    end
    F = fit(  transpose( x_array) , plot_array, 'exp1');

%     semilogy(x_array,plot_array)
%     hold on
%     semilogy(x_array,F(x_array) )
%     hold off

    C =  coeffvalues(F);
    A = C(1);
    xi_inv = - C(2);

     
 end