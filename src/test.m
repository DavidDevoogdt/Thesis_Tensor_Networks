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


M = 10;
N = 10
beta = 0.1;
truncdim = 100;

mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
mpo_base_01 = mpo_base.type_01(0);

mpo_01_M = mpoProduct(mpo_base_01,M,N,100,1);
mpo_02_M = mpoProduct(mpo_base_01,M,N,110,1);

err = error_eigenvalue(mpo_01_M, mpo_02_M,d)

