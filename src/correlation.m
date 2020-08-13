


a=0.01;

%setup T range
T_num = 10; 
T_min =10;
T_max = 100;

T_index = 1:T_num;
T_values = T_min:  (T_max-T_min)/(T_num-1)  :T_max;



%setup m range
m_num = 2; 
m_min =0.01;
m_max = 0.02;

m_index = 1:m_num;
m_values = m_min:  (m_max-m_min)/(m_num-1)  :m_max;


%remembers all values: 
plotstructure = zeros(m_num,T_num,2); 



for T_i = T_index
   for m_i = m_index
        plotstructure(m_i,T_i,:) = fit_correlation_mpo(m_i,T_i,a);
   end
end

for m_i=m_index

    f1 = figure(2*m_i);
    plot( T_values, plotstructure(m_i,:,2))
    hold 
    plot( T_values, T_values*pi/4)
    title('1/xi')
    legend('simulation','analytical')
    hold off

    f2 = figure(2*m_i+1);

    plot( T_values, plotstructure(m_i,:,1))
    hold 
    plot( T_values, 0.8587*(T_values./J).^(1/4))
    legend('simulation','analytical')
    title('A')
    hold off

end

%%
function [A,inv_xi] = fit_correlation_mpo(m,T,a)
    d = 2; % d
    % pre setup
    S_x = [0,1;1,0];
    S_y = [0,-1i;1i,0];
    S_z = [1,0;0,-1];
    I_tensor = eye(2);
    %see http://qpt.physics.harvard.edu/c14.pdf for defintion parameters and
    %hamiltonian

    J=1/2*a;
    %calc properties
    g_c = 1;
    g = g_c-m/(2*J);

    %constuction H tensor final numbering legs: (i1,i2, ..j1,j2)
    H_tensor = -J*( g/2*( ncon( {S_x,I_tensor}, {[-1,-3],[-2,-4]})... 
                   + ncon( {I_tensor,S_x}, {[-1,-3],[-2,-4]}) )...
                   + ncon( {S_z,S_z}, {[-1,-3],[-2,-4]}) );

    H_tensor = complex(H_tensor);

    mpo_gen = generateMPO(d,-H_tensor/T);
    O = mpo_gen.type_01(1);
    % todo make efficient

    N=4;

    x_array = 1:N;
    plot_array = zeros(N,1);
    for i = x_array
        plot_array(i) =  correlation(O,i,S_z);
    end
    F = fit(transpose( x_array) , plot_array, 'exp1');

    %semilogy(plot_array)
    %hold on
    %semilogy(x_array,F(x_array) )
    %hold off


    C =  coeffvalues(F);
    A = C(1);
    inv_xi = - C(2);

end

 %% 
 function c = correlation(O,r,M,A)
     %calculates <psi_0|M_i M_(i+r) exp(beta H)|psi_0>
     %with A MPS corresponding to psi_0
     
     %todo implement this
    
    

     %O_n = Traced_O^(r-1);
     %c =1/(rho_left*rho_right*lambda^(r+1))* ncon( {rho_left,O,M,O_n,O,M,rho_right}, {  [-1,1], [1, 2,3,4],[3,2],[4,5], [5,6,7,8],[7,6],[8]} );
 
 end
 
