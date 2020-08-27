function test
    mpo_type_comparison_exact_generic_order;
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
    
    
    %change these parmas
    
    Order_arr = [1,2,3,4,5,6];% 2:2:6;% [4,6,8]; %keep lower than 8
    order_size = size(Order_arr,2);
    line_spec = ["-"    ,"--"   ,"-."   ,":"    ,"-"    ,"--"   ,"-."   ,":"];
    alphas= [1     ,1      ,1      ,1      ,0.5    ,0.5    ,0.5    ,0.5];
    colors = {[1 0 0 0],[0 1 0 0],[0 0 1 0]};
    legend_Arr = cell(3*order_size,1);
    legend_Arr(:) = {"todo"};
    M=11;
    
    beta_arr = 10.^(  -5:0.1:log10(20) );     
    beta_len = size(beta_arr,2);
    
    MPO_beta_N = cell(beta_len,1);
    for i=1:beta_len
        beta = beta_arr(i);
        MPO_beta_N{i} = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor );
    end
    
    
    
    %
    compare_opts.ref=2;
  
    options5.method="diag";
    options5.testing=0;
    
    opt4.method = "svd";
    opt4.testing = 0;
    opt4.testN=3;
    
    %plot stuff
    plot_counter = 1;
    hold off
    f1=figure(1);
   
    
    for j = 1:order_size
        Order = Order_arr(j);
        plot_structure = cell( 5,beta_len);

        for i = 1:beta_len
            beta = beta_arr(i);
            mpo_base = MPO_beta_N{i};
            
            mpo_base_matrix = mpo_base.H_exp(M-1,1); %buffered, not calculted again every time

            try
                [N5,mpo_N_05] = mpo_base.type_05(Order,options5);
                err_05 = error_eigenvalue({N5,mpo_N_05},"O",mpo_base_matrix,"array",M,d,compare_opts);
                plot_structure{5,i} = err_05; 
            catch
                 err_05=Inf;
                 plot_structure{5,i} = err_05;
            end
            
            try
                [N4,mpo_N_04] = mpo_base.type_04(Order,opt4);
                err_04 = error_eigenvalue({N4,mpo_N_04},"O",mpo_base_matrix,"array",M,d,compare_opts);
                plot_structure{4,i} = err_04;
            catch
                 err_04=Inf;
                 plot_structure{4,i} = err_04;
            end
            
            if Order < 8
                try
                    [N3,mpo_N_03] = mpo_base.type_03(Order,0);
                    err_03 = error_eigenvalue({N3,mpo_N_03},"O",mpo_base_matrix,"array",M,d,compare_opts);
                    plot_structure{3,i} = err_03;
                catch
                    err_03 =Inf;
                    plot_structure{3,i} = err_03;
                end
            end
           

            fprintf("M %d beta %.4e order %d err_03 %.4e err_04 %.4e %.4e \n",M,beta,Order,abs(err_03),abs(err_04),abs(err_05) );

        end

        
       
        
        if Order < 8
            colour = colors{1};
            colour(4) = alphas(j);
            loglog( beta_arr, abs(cell2mat(plot_structure(3,:))), "LineStyle", line_spec(j),"Color",colour );
            legend_Arr{plot_counter}= sprintf("type 03 Order %d",Order );
            plot_counter = plot_counter+1;
            hold on
        end
       
        colour = colors{2};
        colour(4) = alphas(j);
        loglog( beta_arr, abs(cell2mat(plot_structure(5,:))),"LineStyle", line_spec(j),"Color",colour    );
        legend_Arr{plot_counter}= sprintf("type 05 Order %d",Order );
        plot_counter = plot_counter+1;
        
        
        colour = colors{3};
        colour(4) = alphas(j);
        loglog( beta_arr, abs(cell2mat(plot_structure(4,:))),"LineStyle", line_spec(j),"Color",colour );
        legend_Arr{plot_counter}= sprintf("type 01 Order %d",Order );
        plot_counter = plot_counter+1;
        
        title( sprintf("traverse ising g=%.2f",g)   )
        xlabel('$  \beta \cdot J$','Interpreter','latex')
        ylabel("relative err")
        legend(legend_Arr,'Location','northwest')
  
        figure(gcf)
        
        

    end  


    hold off
             
end

