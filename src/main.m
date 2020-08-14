a=0.1;
m=1;
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
H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});
H_1_tensor = -J*g*S_x;


M = 6;

for beta = 0.5:0.5:2
    mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
    mpo_base_01 = mpo_base.type_01(0);
    mpo_base_02 = mpo_base.type_02(0);
    
    
    for N = 2:3
        mpo_N = generateMPO(d,-N*beta*H_1_tensor,-N*beta*H_2_tensor);
        
        mpo_N_01 = mpo_N.type_01(0);
        mpo_N_02 = mpo_N.type_02(0);
        
        mpo_base_N_01 = mpoProduct(mpo_base_01,N);
        mpo_base_N_02 = mpoProduct(mpo_base_02,N);
        
        err_01 = error_eigenvalue(mpo_N_01, mpo_base_N_01,M,d);
        err_02 = error_eigenvalue(mpo_N_02, mpo_base_N_02,M,d);
        
        fprintf("M %d beta %.02f N %d err01 %.4e, err_02 %.4e \n",M,beta,N,abs(err_01), abs(err_02))
    end
end


