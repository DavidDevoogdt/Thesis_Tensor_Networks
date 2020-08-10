
d = 2; % d


delta = 1/2; % from xxz model

timestep = 1;

%constuction H tensor final numbering legs: (i1,i2, ..j1,j2)

% basis chosen such that diagonal form persist for diagonalisation
Hij = [delta/4   0           0           0 
       0         -delta/4    1/2         0
       0         1/2         -delta/4    0
       0         0           0           delta/4].*timestep; 
   
H_tensor = reshape(Hij,  [d,d,d,d]);
I_tensor = eye(d);

mpo_gen = generateMPO(d, H_tensor);
O = mpo_gen.type_01(1);
