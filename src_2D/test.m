function test
    test_H_exp
    test_PEPO
end

function test_H_exp
    d=2;
    S_x = 0.5* [0,1;1,0];
    S_y = 0.5* [0,-1i;1i,0];
    S_z = 0.5* [1,0;0,-1];
    I_tensor = eye(2);

    %ij indices in front, others dont matter hare
    J=1;
    g=2.1;
    
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
   
end

function test_PEPO
    d=2;
    
    %hamiltonian setup
    S_x = 0.5* [0,1;1,0];
    S_y = 0.5* [0,-1i;1i,0];
    S_z = 0.5* [1,0;0,-1];
    I_tensor = eye(2);

    J=1;
    g=2.1;
    
    H_1_tensor = -J*g*reshape( S_x, [d,d]);
    H_2_tensor = -J* reshape( ncon( {S_z,S_z}, {[-1,-3],[-2,-4]}), [d,d,d,d]); 

    opts.testing=0;
    opts.visualise=0;
               
    pos_map = [1,0,0,1;
               1,1,1,1;
               1,1,1,1];

              
    
    map = PEPO.create_map(pos_map);
    map_arg = struct("map", map  );     
    beta_arr = 10.^(-5:0.2:2);
    beta_len = size(beta_arr,2);
    err_arr = zeros( beta_len ,1); 

    for i=1:beta_len
        beta = beta_arr(i);
        pepo = PEPO(d,beta*H_1_tensor,beta*H_2_tensor,2,1,opts);
        [err,prefact] = pepo.calculate_error(map_arg );
        fprintf(" beta %.4e rel err %.4e abs err %.4e \n",beta,abs(err), abs(err)* prefact );
        err_arr(i) = abs(err);

    end
    
    figure();
    loglog(  beta_arr,err_arr );
    title(  sprintf("2D transverse ising M=%d",map.N));
    xlabel( "$\beta \cdot J$","Interpreter","Latex");
    ylabel( "relative error");
    

end