


a=0.1;

%setup T range
T_num = 10; 
T_min = 10;
T_max = 20;

T_index = 1:T_num;
T_values = T_min:  (T_max-T_min)/(T_num-1)  :T_max;



%setup m range
m_num = 2; 
m_min = -0.2;
m_max = -0.1;

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
    title(sprintf('1/xi %d',m_i))
    legend('simulation','analytical')

    hold off

    f2 = figure(2*m_i+1);

    plot( T_values, plotstructure(m_i,:,1))
    hold 
    plot( T_values, 0.8587*(T_values./J).^(1/4))
    legend('simulation','analytical')
    title(sprintf('A %d',m_i))
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

    J=1/(2*a);
    %calc properties
    g_c = 1;
    g = g_c-m/(2*J);

    %constuction H tensor final numbering legs: (i1,i2, ..j1,j2)
    H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});
    H_1_tensor = -J*g*S_x;
               


    mpo_gen = generateMPO(d,-1/T*H_1_tensor,-1/T*H_2_tensor);
    O = mpo_gen.type_01(1);
    
    %todo calc ground state mps for mpo

    N=5;

    x_array = 1:N;
    plot_array = zeros(N,1);
    for i = x_array
        plot_array(i)=  thermal_correlation(O,i,S_z);
    end
    F = fit(  transpose( x_array) , plot_array, 'exp1');

    semilogy(x_array,plot_array)
    hold on
    semilogy(x_array,F(x_array) )
    hold off


    C =  coeffvalues(F);
    A = C(1);
    inv_xi = - C(2);

end

 %% 
 function c = thermal_correlation(O,r,M)

    
     
    bond_d  = size(O,1);
    d =size(O,2);
    
    I =eye(d);
    
    Traced_O = ncon( {O},{[-1,1,1,-2]});
    O_n = Traced_O^(r-1);
    
%     [V,D,W] = eig(Traced_O);
%     rho_right= V(:,1);
%     rho_left = W(:,1)';
%     lambda = D(1,1);
%     O_n = Traced_O^(r-1);
%     c =1/(rho_left*rho_right*lambda^(r+1))* ncon( {rho_left,O,M,O_n,O,M,rho_right}, {  [-1,1], [1, 2,3,4],[3,2],[4,5], [5,6,7,8],[7,6],[8,-2]} );
%  
    end_vector = zeros( bond_d,1);
    end_vector(1) = 1;
    
    
    
    c =ncon( {end_vector,O,M,O_n,O,M,end_vector}, {  [1], [1, 2,3,4],[3,2],[4,5], [5,6,7,8],[7,6],[8,]} )/...
       ncon( {end_vector,O,I,O_n,O,I,end_vector}, {  [1], [1, 2,3,4],[3,2],[4,5], [5,6,7,8],[7,6],[8]} );
 
     
 end
 
