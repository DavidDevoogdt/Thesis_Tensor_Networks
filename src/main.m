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
H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});    %2 site operator
H_1_tensor = -J*g*S_x;                                  %on every sing site


M = 100;
beta = 0.1;
truncdim = 100;

mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
mpo_base_01_M = mpoProduct(mpo_base.type_01(0),M,1,100);
mpo_base_02_M = mpoProduct(mpo_base.type_02(0),M,1,400);

for N = [5,10,15]
    
    mpo_N = generateMPO(d,-beta/N*H_1_tensor,-beta/N*H_2_tensor );

    mpo_N_01 = mpo_N.type_01(0);
    mpo_N_02 = mpo_N.type_02(0);
    
    mpo_01_M = mpoProduct(mpo_N_01,M,N,truncdim);
    err_01 = error_eigenvalue(mpo_base_01_M, mpo_01_M,d);
    
    mpo_02_M = mpoProduct(mpo_N_02,M,N,truncdim);
    err_02 = error_eigenvalue(mpo_base_02_M, mpo_02_M,d);

    err = error_eigenvalue(mpo_01_M, mpo_02_M,d);
    
    fprintf("truncdim %3d M %d beta %.02f N %3d err_01 %.4e err_02 %.4e err %.4e\n",truncdim,M,beta,N,abs(err_01),abs(err_02),abs(err));
end

% for beta = 1:3
%     mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
%     mpo_base_01 = mpo_base.type_01(0);
%     mpo_base_02 = mpo_base.type_02(0);
% 
%     for N = 10:10
%         mpo_N = generateMPO(d,-beta/N*H_1_tensor,-beta/N*H_2_tensor);
%         
%         mpo_N_01 = mpo_N.type_01(0);
%         mpo_N_02 = mpo_N.type_02(0);
%         
%         mpo_N_01_M = mpoProduct(mpo_N_01,M,N,100);
%         mpo_N_02_M = mpoProduct(mpo_N_02,M,N,100);
%         
%         mpo_base_01_M = mpoProduct(mpo_base_01,M,1,200);
%         mpo_base_02_M = mpoProduct(mpo_base_02,M,1,200);
%         
%         err_01 = error_eigenvalue(mpo_N_01_M, mpo_base_01_M,d);
%         err_02 = error_eigenvalue(mpo_N_02_M, mpo_base_02_M,d);
%         
%         fprintf("M %d beta %.02f N %d err01 %.4e, err_02 %.4e \n",M,beta,N,abs(err_01), abs(err_02))
%     end
% end


