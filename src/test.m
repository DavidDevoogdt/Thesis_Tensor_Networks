function test
    %test_truncateO();
    %test_mpo_prod_2
    mpo_type_comparison();
end

%90 is sufficient
function test_truncateO
    a=1;
    m=0.7;
    d = 2; 
    
    S_x = [0,1;1,0];
    S_y = [0,-1i;1i,0];
    S_z = [1,0;0,-1];
    I_tensor = eye(2);
    
    J=1/(2*a);
    %calc properties
    g_c = 1;
    g = g_c-m/(2*J);
    
    beta=1.0;
    truncdim=90; %seems about minimum for mpo type 2
    M=10;

    %constuction H tensor final numbering legs: (i1,i2, ..j1,j2)
    H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});    %2 site operator
    H_1_tensor = -J*g*S_x;                                  %on every sing site

    mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
    O = mpo_base.type_02(0);
    

    mpo_01 = generateList(O,M,truncdim,1); %symmetrically truncated
    mpo_02 = generateList(O,M,-1); %untruncated list
    
    err = error_eigenvalue(mpo_01, mpo_02,d);
    fprintf("truncdim %d err = %.4e\n",truncdim,abs(err));
    
end

function test_mpo_prod_2
    a=1;
    m=0.7;
    d = 2; 
    
    S_x = [0,1;1,0];
    S_y = [0,-1i;1i,0];
    S_z = [1,0;0,-1];
    I_tensor = eye(2);
    
    J=1/(2*a);
    %calc properties
    g_c = 1;
    g = g_c-m/(2*J);
    
    beta=1.0;
    truncdim=150; %seems about minimum for mpo type 2
    M=10;

    %constuction H tensor final numbering legs: (i1,i2, ..j1,j2)
    H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});    %2 site operator
    H_1_tensor = -J*g*S_x;                                  %on every sing site

    mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
    O = mpo_base.type_01(0);
    list = generateList(O,M,-1);

    %mpoProduct(O_cell,N,truncdim,testing,type)
    
    mpo_01 = mpoProduct(list,4,truncdim,1,2); %
    mpo_02 = mpoProduct(list,4,truncdim+10,1,2);%
    
    err = error_eigenvalue(mpo_01, mpo_02,d);
    fprintf("truncdim %d err = %.4e\n",truncdim,abs(err));
    
end

function mpo_type_comparison
    a=1;
    m=0.4;
    d = 2; % d
    % pre setup
    S_x = [0,1;1,0];
    S_y = [0,-1i;1i,0];
    S_z = [1,0;0,-1];
    I_tensor = eye(2);

    J=1;
    g=0.6;


   
    H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});    %2 site operator
    H_1_tensor = -J*g*S_x;                                  %on every sing site


    M = 10;
    beta = 0.1;
    truncdim = 200;



    mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
    mpo_base_01_M = generateList(mpo_base.type_01(0),M,-1); %not truncated
    mpo_base_02_M = generateList(mpo_base.type_02(1),M,90); %truncdim 90 does not make mpo worse

    for N = [2,3]
        mpo_N = generateMPO(d,-beta/N*H_1_tensor,-beta/N*H_2_tensor );


        fprintf('N=%3d started at %s \n', N, datestr(now,'HH:MM:SS.FFF'))

        mpo_N_01_cell = truncateO(mpo_N.type_01(0),-1);
        mpo_N_01 = mpoProduct(mpo_N_01_cell,N,truncdim,1);
        mpo_N_01_M = generateList(mpo_N_01,M,-1);

        err_01 = error_eigenvalue(mpo_base_01_M, mpo_N_01_M ,d);
        fprintf(' type 01 finished at %s \n',datestr(now,'HH:MM:SS.FFF'))

        
        fprintf("truncdim %3d M %d beta %.02f N %3d err_01 %.4e \n",truncdim,M,beta,N,abs(err_01));
        
%         mpo_N_02_cell = truncateO(mpo_N.type_02(0),90);
%         mpo_N_02 = mpoProduct(mpo_N_02_cell,N,truncdim,1);
%         mpo_N_02_M = generateList(mpo_N_02,M,-1);
% 
%         err_02 = error_eigenvalue(mpo_base_02_M, mpo_N_02_M ,d);
%         fprintf(' type 02 finished at %s \n',datestr(now,'HH:MM:SS.FFF'))
% 
% 
%         err = error_eigenvalue(mpo_N_01_M, mpo_N_02_M,d);
% 
%         fprintf("truncdim %3d M %d beta %.02f N %3d err_01 %.4e err_02 %.4e err %.4e\n",truncdim,M,beta,N,abs(err_01),abs(err_02),abs(err));
    end

end

