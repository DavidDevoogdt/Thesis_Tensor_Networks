function test
    
    %test_trunc ateO();
    %test_mpo_prod_2
    %mpo_type_comparison();
    %mpo_type_comparison_asym
    %compare_with_exact_hamiltonian
    %compare_01_with_exact_hamiltonian_assym
    %%compare_01_with_exact_hamiltonian_sym
    %mpo_type_comparison_exact
    %mpo_type_comparison_exact_2
    mpo_type_comparison_exact_generic_order
    %test_type_02_decomp
    %H_exp_precision_test
end


%makes it worse
function H_exp_precision_test
    a=1;
    m=0.4;
    d = 2; % d
    % pre setup
    S_x = 0.5* [0,1;1,0];
    S_y = 0.5* [0,-1i;1i,0];
    S_z = 0.5* [1,0;0,-1];
    I_tensor = eye(2);


    J=1;
    g=0.7;
   
    H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});    %2 site operator
    H_1_tensor = -J*g*S_x;                                  %on every sing site

    M_arr = 6:11;
    
    m_size = size(M_arr,2);


    
    for j = 1:m_size
        M = M_arr(j);

        beta_arr = [5.0,10,20,50];
       
        
        len = size(beta_arr,2);
        
        N_array = 1:5;
        n_len =  size(beta_arr,2);
        
        mpo_array = cell(n_len,1);

        for i = 1:len
            beta = beta_arr(i);

            
            
            for n_index=1:n_len
                N = N_array(n_index);
                
                mpo_N = generateMPO(d,-beta/N*H_1_tensor,-beta/N*H_2_tensor );
                mpo_N_matrix = mpo_N.H_exp(M-1,1);
                
                prefact = trace(mpo_N_matrix);
                mpo_matrix = prefact^N.*  (mpo_N_matrix./prefact)^N;
                
                mpo_array{n_index} = mpo_matrix;
                %fprintf("N %d beta %f err %.4e\n",N,beta, eigs(Z,1));
                
            end
            
            for N = 1:n_len-1
            
                fprintf("beta %f err H_%d H_%d %.4e\n",beta,N,N+1, eigs( mpo_array{N}-  mpo_array{N+1} ,1));
            end
        end
 
    end
  

end


function test_type_02_decomp
     a=1;
    m=0.4;
    d = 2; % d
    % pre setup
    S_x = [0,1;1,0];
    S_y = [0,-1i;1i,0];
    S_z = [1,0;0,-1];
    I_tensor = eye(2);

    J=1;
    delta=0.5;

    %H_2_tensor = ncon( {S_x,S_x}, {[-1,-3],[-2,-4]})...
    %            +ncon( {S_y,S_y}, {[-1,-3],[-2,-4]})...
    %      +delta*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});%2 site operator
    %H_1_tensor = 0;                                  %on every sing site

     
    J=1;
    g=0.6;
   
    H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});    %2 site operator
    H_1_tensor = -J*g*S_x;                                  %on every sing site

        
    %beta_arr = [0.1,0.5,1.0,2.0,4.0,8.0];
    beta_arr =[2];

    len = size(beta_arr,2);



    for i = 1:len
        beta = beta_arr(i);
        mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
        
        fprintf("type_02\n");
        mpo_base.type_02(1);
        
        fprintf("type_02\n");
        mpo_base.type_03(1);%just interested in printed decomp errors

    end

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
    O = mpo_base.type_02(1);
    

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

function mpo_type_comparison_asym
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


    M = 8;
    beta = 0.1;
    truncdim = 100;

    mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
    mpo_base_01_M = generateList(mpo_base.type_01(0),M,-1); %not truncated
    mpo_base_02_M = generateList(mpo_base.type_02(0),M,90); %truncdim 90 does not make mpo worse

    for N = [2]
        mpo_N = generateMPO(d,-beta/N*H_1_tensor,-beta/N*H_2_tensor );


        fprintf('N=%3d started at %s \n', N, datestr(now,'HH:MM:SS.FFF'))

        mpo_01_M = generateList( mpo_N.type_01(0),M,-1);
        mpo_N_01_M =  mpoProduct(mpo_01_M,N,truncdim,1,2);

        err_01 = error_eigenvalue(mpo_base_01_M, mpo_N_01_M ,d);
        fprintf(' type 01 finished at %s \n',datestr(now,'HH:MM:SS.FFF'))
        
       

        mpo_02_M = generateList( mpo_N.type_02(0),M,-1);
        mpo_N_02_M =  mpoProduct(mpo_02_M,N,truncdim,1,2);

        err_02 = error_eigenvalue(mpo_base_02_M, mpo_N_02_M ,d);
        fprintf(' type 02 finished at %s \n',datestr(now,'HH:MM:SS.FFF'))

        err = error_eigenvalue(mpo_N_01_M, mpo_N_02_M,d);
 
        fprintf("truncdim %3d M %d beta %.02f N %3d err_01 %.4e err_02 %.4e err %.4e\n",truncdim,M,beta,N,abs(err_01),abs(err_02),abs(err));
 
    end

end

function mpo_type_comparison_exact
    a=1;
    m=0.4;
    d = 2; % d
    % pre setup
    S_x = 0.5* [0,1;1,0];
    S_y = 0.5* [0,-1i;1i,0];
    S_z = 0.5* [1,0;0,-1];
    I_tensor = eye(2);

    J=1;
    delta=0.5;

%     H_2_tensor = ncon( {S_x,S_x}, {[-1,-3],[-2,-4]})...
%              S  +ncon( {S_y,S_y}, {[-1,-3],[-2,-4]})...
%           +delta*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});%2 site operator
%     H_1_tensor = zeros(d,d);                                  %on every sing site

     
    J=1;
    g=0.6;
   
    H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});    %2 site operator
    H_1_tensor = -J*g*S_x;                                  %on every sing site

    M_arr = 6:11;
    
    m_size = size(M_arr,2);
    
    m_struct = cell(m_size,2,2  );

    
    for j = 1:m_size
        M = M_arr(j);
        
        %beta_arr = 4:8;
        %beta_arr = [0.1,0.5,1,2,4,6,7,7.5,8.1,10];
        %fit_inidces = 5:10;
        %beta_arr = [0.01,0.05,0.1,0.2, 0.5,1];
        %fit_inidces = 1:6;

        %beta_arr = [0.001,0.002,0.005,0.01,0.02,05,0.1,0.2, 0.5,1.0];
        beta_arr = 2:2:8
        fit_inidces = 1:6;
        
        len = size(beta_arr,2);

        plot_structure = cell( 3,len);

        for i = 1:len
            beta = beta_arr(i);

            mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
            mpo_base_matrix = mpo_base.H_exp(M-1,1);

            mpo_N_01_M = generateList(mpo_base.type_01(0),M,-1);
            err_01 = error_eigenvalue(mpo_N_01_M,mpo_base_matrix ,d);

            mpo_N_02_M = generateList(mpo_base.type_02(0),M,-1);
            err_02 = error_eigenvalue( mpo_N_02_M,mpo_base_matrix,d);
            
            mpo_N_03_M = generateList(mpo_base.type_03(0),M,-1);
            err_03 = error_eigenvalue( mpo_N_03_M,mpo_base_matrix,d);
            
            err_04 = error_eigenvalue( mpo_N_02_M,mpo_N_03_M,d);

            plot_structure{1,i} = err_01;
            plot_structure{2,i} = err_02;
            plot_structure{3,i} = err_03;

            fprintf("M %d beta %.02f err_01 %.4e err_02 %.4e err_03 %.4e mpo02-mpo03 err %.4e \n",M,beta,abs(err_01),abs(err_02),abs(err_03),abs(err_04) );

        end

        
        y=(abs(  cell2mat(plot_structure(1,:))));

        
        F1 = fit( beta_arr(fit_inidces).' ,y(fit_inidces).' , 'exp1');
        C1 =  coeffvalues(F1);
        A1 = C1(1);
        n1 = C1(2);

        y=(abs(  cell2mat(plot_structure(2,:))));
        
        F2 = fit(  beta_arr(fit_inidces).' , y(fit_inidces).' , 'exp1');
        C2 =  coeffvalues(F2);
        A2 = C2(1);
        n2 = C2(2);
        
        y=(abs(  cell2mat(plot_structure(3,:))));
        
        F3 = fit(  beta_arr(fit_inidces).' , y(fit_inidces).' , 'exp1');
        C3 =  coeffvalues(F3);
        A3 = C3(1);
        n3 = C3(2);

        figure(1)
        semilogy( beta_arr,abs(  cell2mat(plot_structure(1,:) )), beta_arr, F1(beta_arr),beta_arr, abs(cell2mat(plot_structure(2,:))),beta_arr, F2(beta_arr),beta_arr, abs(cell2mat(plot_structure(3,:))),beta_arr, F3(beta_arr) );
        legend("type 01","01 fit","type 02","02 fit","type 03","03 fit",'Location','northwest')
        title(M)
        xlabel("beta")
        ylabel("err")

   
        
        m_struct{j,1,1} = A1;
        m_struct{j,1,2} = A2;
        m_struct{j,1,3} = A3;
        
        m_struct{j,2,1} = n1;
        m_struct{j,2,2} = n2;
        m_struct{j,2,3} = n3;
    end
  
    figure(1)
    semilogy(M_arr,cell2mat(m_struct(:,1,1)), M_arr,cell2mat(m_struct(:,1,2)),M_arr,cell2mat(m_struct(:,1,3)) );
    legend("type 01","type 02","type 03");
    xlabel("M")
    ylabel("A")
    
    title("$ error= A \cdot \beta^{ n} $ ",'Interpreter','latex')
    
        
    figure(2)
    plot( M_arr,cell2mat(m_struct(:,2,1)), M_arr,cell2mat(m_struct(:,2,2)), M_arr,cell2mat(m_struct(:,2,3)) );
    legend("type 01","type 02","type 03");
    xlabel("M")
    ylabel("n")
    
    title("$ error= A \cdot \beta^{ n} $ ",'Interpreter','latex')
    
    
    
end

function mpo_type_comparison_exact_2
    a=1;
    m=0.4;
    d = 2; % d
    % pre setup
    S_x = 0.5* [0,1;1,0];
    S_y = 0.5* [0,-1i;1i,0];
    S_z = 0.5* [1,0;0,-1];
    I_tensor = eye(2);

    J=1;
    delta=0.5;

%     H_2_tensor = ncon( {S_x,S_x}, {[-1,-3],[-2,-4]})...
%              S  +ncon( {S_y,S_y}, {[-1,-3],[-2,-4]})...
%           +delta*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});%2 site operator
%     H_1_tensor = zeros(d,d);                                  %on every sing site

     
    J=1;
    g=2.1;
   
    H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});    %2 site operator
    H_1_tensor = -J*g*S_x;                                  %on every sing site

    M_arr = 6:11;
    
    m_size = size(M_arr,2);
    

    opts.ref=2;
  

    
    hold off

    f1=figure(1);
    

    
    for j = 1:m_size
        M = M_arr(j);
        
        %beta_arr = 4:8;
        %beta_arr = [0.1,0.5,1,2,4,6,7,7.5,8.1,10];
        %fit_inidces = 5:10;
        %beta_arr = [0.01,0.05,0.1,0.2, 0.5,1];
        %fit_inidces = 1:6;

        %beta_arr = [0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2, 0.5,1.0,2.0,5.0,10,20,50];
        beta_arr = 10.^(  -3:0.05:1.5  );

        
        len = size(beta_arr,2);

        plot_structure = cell( 4,len);

        for i = 1:len
            beta = beta_arr(i);

            mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
            mpo_base_matrix = mpo_base.H_exp(M-1,1);

            [N1,mpo_N_01] = mpo_base.type_01(0);
            err_01 = error_eigenvalue( {N1,mpo_N_01},"O",mpo_base_matrix,"array",M,d,opts);

            [N2,mpo_N_02] = mpo_base.type_02(0);
            err_02 = error_eigenvalue({N2,mpo_N_02},"O",mpo_base_matrix,"array",M,d,opts);
            
            [N3,mpo_N_03] = mpo_base.type_03(0);
            err_03 = error_eigenvalue({N3,mpo_N_03},"O",mpo_base_matrix,"array",M,d,opts);
            
            [N4,mpo_N_04] = mpo_base.type_04(6,1);
            err_04 = error_eigenvalue({N4,mpo_N_04},"O",mpo_base_matrix,"array",M,d,opts);
            
            
            plot_structure{1,i} = err_01;
            plot_structure{2,i} = err_02;
            plot_structure{3,i} = err_03;
            plot_structure{4,i} = err_04;


            fprintf("M %d beta %.4e err_01 %.4e err_02 %.4e err_03 %.4e err_04 %.4e \n",M,beta,abs(err_01),abs(err_02),abs(err_03),abs(err_04) );

        end

        
        loglog( beta_arr,abs(  cell2mat(plot_structure(1,:) )), "color","blue");
        hold on
        loglog( beta_arr, abs(cell2mat(plot_structure(2,:))),"color","red");
        
        loglog( beta_arr, abs(cell2mat(plot_structure(3,:))),"color","green");
        loglog( beta_arr, abs(cell2mat(plot_structure(4,:))),"color","cyan");
        
        
        legend("type 01","type 02","type 03","type 04",'Location','northwest')
         xlabel("beta")
         ylabel("err")

  
        figure(gcf)
        
        

    end  
    hold off
             
end



function mpo_type_comparison_exact_generic_order
    a=1;
    m=0.4;
    d = 2; % d
    % pre setup
    S_x = 0.5* [0,1;1,0];
    S_y = 0.5* [0,-1i;1i,0];
    S_z = 0.5* [1,0;0,-1];
    I_tensor = eye(2);

    J=1;
    delta=0.5;

%     H_2_tensor = ncon( {S_x,S_x}, {[-1,-3],[-2,-4]})...
%              S  +ncon( {S_y,S_y}, {[-1,-3],[-2,-4]})...
%           +delta*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});%2 site operator
%     H_1_tensor = zeros(d,d);                                  %on every sing site

     
    J=1;
    g=2.1;
   
    H_2_tensor = -J*ncon( {S_z,S_z}, {[-1,-3],[-2,-4]});    %2 site operator
    H_1_tensor = -J*g*S_x;                                  %on every sing site

    Order_arr = [2,4,6,8];
    
    order_size = size(Order_arr,2);
    

    opts.ref=2;
  
    legend_Arr = cell(order_size*2,1);
    legend_Arr(:) = {"todo"};
    
    hold off

    f1=figure(1);
   
    

    
    for j = 1:order_size
        Order = Order_arr(j);
        
        M=10;
        
      
        %beta_arr = [0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2, 0.5,1.0,2.0,5.0,10,20,50];
        beta_arr = 10.^(  -3:0.2:1.5  );

        
        len = size(beta_arr,2);

        plot_structure = cell( 4,len);

        for i = 1:len
            beta = beta_arr(i);

            mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
            mpo_base_matrix = mpo_base.H_exp(M-1,1);

            [N3,mpo_N_03] = mpo_base.type_03(0);
            err_03 = error_eigenvalue({N3,mpo_N_03},"O",mpo_base_matrix,"array",M,d,opts);
            
            [N4,mpo_N_04] = mpo_base.type_04(Order,1);
            err_04 = error_eigenvalue({N4,mpo_N_04},"O",mpo_base_matrix,"array",M,d,opts);
            
            
            %plot_structure{1,i} = err_01;
            %plot_structure{2,i} = err_02;
            plot_structure{3,i} = err_03;
            plot_structure{4,i} = err_04;


            fprintf("M %d beta %.4e order %d err_03 %.4e err_04 %.4e \n",M,beta,Order,abs(err_03),abs(err_04) );

        end

        
        %loglog( beta_arr,abs(  cell2mat(plot_structure(1,:) )), "color","blue");
        
        %loglog( beta_arr, abs(cell2mat(plot_structure(2,:))),"color","red");
        
        loglog( beta_arr, abs(cell2mat(plot_structure(3,:)))  );
        hold on
        loglog( beta_arr, abs(cell2mat(plot_structure(4,:))));
        
        legend_Arr{2*j-1}= sprintf("type 03" );
        legend_Arr{2*j}= sprintf("type 01 Order %d",Order );
        
         xlabel("beta")
         ylabel("err")
        legend(legend_Arr,'Location','northwest')
  
        figure(gcf)
        
        

    end  


    hold off
             
end



function check_convergence_between_mpo
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
    truncdim = 150;
    beta = 1;
    
    mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
    exact_diag = mpo_base.H_exp(M-1,1);

end

function compare_01_with_exact_hamiltonian_assym
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
    truncdim = 150;
    beta = 1;
    
    mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
    exact_diag = mpo_base.H_exp(M-1,1);

    
    for N=[1,2,3]

        mpo_N = generateMPO(d,-beta/N*H_1_tensor,-beta/N*H_2_tensor );
        mpo_N_M = mpoProduct( generateList(mpo_N.type_01(0),M,-1),N,truncdim,0,2) ; %not truncated
        err_01 = error_eigenvalue(mpo_N_M, exact_diag,d);

       
        fprintf("asym vs diag: truncdim %2d M %d beta %.02f N %2d error %.4e \n",truncdim,M,beta,N,abs(err_01));
    end
end

function compare_01_with_exact_hamiltonian_sym
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
    truncdim = 150;
    beta = 1;
    
    mpo_base = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
    exact_diag = mpo_base.H_exp(M-1,1);

    
    for N=[1,2,3]

        mpo_N = generateMPO(d,-beta/N*H_1_tensor,-beta/N*H_2_tensor );
        mpo = mpoProduct( truncateO(mpo_N.type_01(0),-1) ,N,truncdim,0) ; %not truncated
        err_01 = error_eigenvalue( generateList(mpo,M,-1), exact_diag,d);
       
        fprintf(" sym vs diag: truncdim %2d M %d beta %.02f N %2d error %.4e \n",truncdim,M,beta,N,abs(err_01));
    end
end

