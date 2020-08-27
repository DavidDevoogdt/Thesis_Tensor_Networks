function test
test_H_exp
end

function test_H_exp


    d=2;
    S_x = 0.5* [0,1;1,0];
    S_y = 0.5* [0,-1i;1i,0];
    S_z = 0.5* [1,0;0,-1];
    I_tensor = eye(2);


    %ij indices in front, others dont matter hare
    H_1_tensor =  reshape(0.3* S_x, [d,d,1,1,1,1]);
    H_2_tensor =  reshape( ncon( {S_z,S_z}, {[-1,-3],[-2,-4]}), [d,d,d,d,1,1,1,1,1,1]); 

    mpo = generatePEPO(d,H_1_tensor,H_2_tensor);

    pos_map = [1,0,0;
               1,0,0;
               1,1,1];

    Z1= mpo.H_exp(pos_map);       


    pos_map = [1,1,0;
               0,1,1;
               0,0,1];

    Z2= mpo.H_exp(pos_map);  
    
    A=Z2-Z1; %should be zero because it the same length chain numbered in same way
end