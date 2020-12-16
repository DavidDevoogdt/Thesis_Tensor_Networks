function test
    %test_H_exp
    test_PEPO
end

function test_H_exp
    fprintf("test_H_exp\n")

    d=2;
    S_x = 0.5* [0,1;1,0];
    S_y = 0.5* [0,-1i;1i,0];
    S_z = 0.5* [1,0;0,-1];
    I_tensor = eye(2);

    %ij indices in front, others dont matter hare
    J=1;
    g=0.4;
    
    opts.testing=0;
    opts.visualise=0;
    
    H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});    %2 site operator
    H_1_tensor = -J*g*S_x;                                  %on every sing site
    
    pepo = PEPO(d,H_1_tensor,H_2_tensor,2,1,opts);

    pos_map = [1,0,0;
               1,0,0;
               1,1,1];

    Z1= pepo.H_exp(struct("pos_map",pos_map));       

    prefact = 1/10;
    Z2= prefact^5* pepo.H_exp(struct("pos_map",pos_map),prefact);  
    
    A=Z2-Z1; %should be close to zero

    pos_map = [1,1,0;
               0,1,1;
               0,0,1];

    Z2=pepo.H_exp(struct("pos_map",pos_map));    
           
    fprintf( "diff= %f \n", sum(  abs( reshape(Z2-Z1,[],1)))) %should be zero because it the same length chain numbered in same way       
   
    fprintf("\n end test_H_exp\n")
    
end

function test_PEPO

    fprintf("\n")

    d=2;
    
    %hamiltonian setup
    S_x =  [0,1;1,0];
    S_y =  [0,-1i;1i,0];
    S_z =  [1,0;0,-1];
    I_tensor = eye(2);

    
    J=1;
    g=0.01;
    
    H_1_tensor = -J*g*S_x;
    H_2_tensor = -J* ( reshape( ncon( {S_z,S_z}, {[-1,-3],[-2,-4]}), [d,d,d,d]));
                        

    opts.testing=0;
    opts.visualise=0;
               
    pos_map = [0,0,1,1;
               1,1,1,1;
               0,0,1,1];


    map = PEPO.create_map(pos_map); 

%     H_1_tensor = 0*eye(d);
%     H_2_tensor = -J* ( reshape( ncon( {S_z,S_z}, {[-1,-3],[-2,-4]}), [d,d,d,d]) +...
%                         0.5*g* reshape( ncon( {S_x,eye(d)}, {[-1,-3],[-2,-4]}), [d,d,d,d])+...
%                         0.5*g* reshape( ncon( {eye(d),S_x}, {[-1,-3],[-2,-4]}), [d,d,d,d]));
%    
    T_c = 2*J/(log(1+sqrt(2)));

    T =2*J./( asinh( (1-(0.7:-0.1:0).^8 ).^(-4))   );
    T=T(T>1);
    
    T = [ 1:0.1:min(T) ,T(T>1), T_c:0.05:(T_c+0.15)]   ;
    
    
    m_arr_theory = T;
    small_T = T<T_c;
    
    m_arr_theory(small_T) = (1-sinh((2*J)./T(small_T)).^(-4)).^(1/8)  ;
    m_arr_theory(~small_T) = 0;

    %beta_arr = 10.^(-2:0.1:0);
    beta_arr=1./T;
    
    beta_len = size(beta_arr,2);
    m_arr = zeros( beta_len ,1); 

    for i=1:beta_len
        beta = beta_arr(i);
        pepo = PEPO(d,-beta*H_1_tensor,-beta*H_2_tensor,3,1,opts);
        %[err,prefact] = pepo.calculate_error(map );
        
        %[A,G,lambda,ctr,error] = pepo.vumps();
        
        mag = pepo.get_expectation( S_z  );
        
                
        fprintf(" T %.4e mag %.4e ,Theory: %.4e \n",1/beta,abs(mag),  m_arr_theory(i));
        m_arr(i) = abs(mag,5);
        
        %fprintf(" beta %.4e rel err %.4e abs err %.4e \n",beta,abs(err), abs(err)* prefact );
        %err_arr(i) = abs(err);

    end
    
    figure();
    %loglog(  beta_arr,err_arr );
    plot(T,m_arr );

    hold on
    
%     
%     title(  sprintf("2D transverse ising M=%d",map.N));
%     xlabel( "$\beta \cdot J$","Interpreter","Latex");
%     ylabel( "relative error");
    
    title(  sprintf("2D transverse Ising, g=%.4f ",g));
    xlabel( "$\frac{k T}{J}$","Interpreter","Latex");
    ylabel( "m");
    

end