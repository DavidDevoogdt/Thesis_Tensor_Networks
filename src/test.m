function test
    mpo_type_comparison_exact_generic_order;
end

function mpo_type_comparison_exact_generic_order

    
    
    
    %change this
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model = "t_ising";
    
    
    simul.Order_arr =  [2,3,4,5,6,7];
    
    simul.types = [3,4];
    simul.M = 8;
    simul.beta_arr = 10.^(  -3:0.5:log10(40) );
    
    generate_opts.testing=0;
    
    opt3 = struct([]);
    
    opt4.method = "svd";
    opt4.to_matrix = 0; %keep in cell form
    
    opt5.method="diag";
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       
       
    S_x_2 = 0.5* [0,1;1,0];
    S_y_2 = 0.5* [0,-1i;1i,0];
    S_z_2 = 0.5* [1,0;0,-1];

    S_x_3 = 1/sqrt(2)* [0,1,0;1,0,1;0,1,0;];
    S_y_3 = 1/sqrt(2)*[0, -1i,0;1i,0, -1i;0, 1i, 0;];
    S_z_3 = [1,0,0;0,0,0;0,0,1;];
       
    switch model
        case "t_ising"
            d=2;
            J=1;
            g=1.05;
            
            H_2_tensor =-J*ncon( {S_z_2,S_z_2}, {[-1,-3],[-2,-4]});
            H_1_tensor =-J*g*S_x_2 ;
    
            opt4.single_threshold = -1;
            opt4.double_threshold = -1;
            
            
            simul.title = sprintf("transverse ising g=%.3f",g);
            
        case "XXZ_3D"
            d = 3; 
            delta = 0.5;

            H_2_tensor =-( ncon( {S_x_3,S_x_3}, {[-1,-3],[-2,-4]})+...
                           ncon( {S_y_3,S_y_3}, {[-1,-3],[-2,-4]})+...
                     delta*ncon( {S_x_3,S_x_3}, {[-1,-3],[-2,-4]}) ...
                          );   
            H_1_tensor =  zeros(d) ;   
            
            simul.title = sprintf("XXZ 3D delta=%.3f",delta);
            
        case "Heisenberg_2D"    
            d=2;
            H_2_tensor = ncon( {S_x_2,S_x_2}, {[-1,-3],[-2,-4]})...
                        +ncon( {S_y_2,S_y_2}, {[-1,-3],[-2,-4]})...
                        +ncon( {S_z_2,S_z_2}, {[-1,-3],[-2,-4]});
            H_1_tensor = zeros(d); 
            
            opt4.single_threshold = 1e-10;
            opt4.double_threshold = 1e-8;
            
            
            
            simul.title = sprintf("Heisenberg");
            
        otherwise
            error("unknown model")
            
    end
    
    %%%plot stuff
    order_size = size(simul.Order_arr,2);
    line_spec = ["-"    ,"--"   ,"-."   ,":"    ,"-"    ,"--"   ,"-."   ,":"];
    alphas= [1     ,1      ,1      ,1      ,0.5    ,0.5    ,0.5    ,0.5];
    colors = {[1 0 0 0],[0 1 0 0],[0 0 1 0]}; %reserved for specific type

    legend_Arr = cell(size(simul.types,2)*order_size,1);
    legend_Arr(:) = {"todo"};
    beta_len = size(simul.beta_arr,2);

    plot_counter = 1;
    hold off
    figure(1);
    %%%
   
    
    %pregenerate all mpo's
    MPO_beta_N = cell(beta_len,1);
    for i=1:beta_len
        beta = simul.beta_arr(i);
        
        %generateMPO(d,H_1_tensor,H_2_tensor,type,order,type_opts,opts)
        MPO_beta_N{i} = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,0 );
    end

    %
    compare_opts.ref=2;
  
    %loop over different orders
    for j = 1:order_size
        Order = simul.Order_arr(j);
        plot_structure = cell( 5,beta_len);
            
        
        
        %loop over temps
        for i = 1:beta_len
            beta = simul.beta_arr(i);
            mpo_base = MPO_beta_N{i};
            mpo_base_matrix = mpo_base.H_exp(simul.M-1,1); %buffered, not calculted again every time

            fprintf("M %d beta %.4e order %d",simul.M,beta,Order);
            
            %loop of simulation types
            for t=1:size(simul.types,2)
               switch simul.types(t)
                   case 3
          
                      mpo_3 = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,3,Order,opt3,generate_opts);
                      err_03 = error_eigenvalue(mpo_3,mpo_base_matrix,"array",simul.M,d,compare_opts);
                      fprintf(" err 03 %.4e",err_03);
                      plot_structure{3,i} = err_03;
                      
                   case 4
                      mpo_4 = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,4,Order,opt4,generate_opts);
                      err_04 = error_eigenvalue(mpo_4,mpo_base_matrix,"array",simul.M,d,compare_opts);
                      fprintf(" err 04 %.4e",err_04);
                      plot_structure{4,i} = err_04;
                   case 5 
                      mpo_5 = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,5,Order,opt5,generate_opts);
                      err_05 = error_eigenvalue(mpo_5,mpo_base_matrix,"array",simul.M,d,compare_opts);
                      fprintf(" err 05 %.4e",err_05);
                      plot_structure{5,i} = err_05; 
                   otherwise
                       error("unknown type")
               end
            end         

            fprintf("\n");

        end
        
        %plotting loop for current order
        for t=1:size(simul.types,2)
               switch simul.types(t)
                   case 3
                      colour = colors{1};
                      colour(4) = alphas(j);
                      loglog( simul.beta_arr, abs(cell2mat(plot_structure(3,:))), "LineStyle", line_spec(j),"Color",colour );
                      legend_Arr{plot_counter}= sprintf("type 03 Order %d",Order );
                      plot_counter = plot_counter+1;
                      hold on
                   case 4
                      colour = colors{3};
                        colour(4) = alphas(j);
                        loglog( simul.beta_arr, abs(cell2mat(plot_structure(4,:))),"LineStyle", line_spec(j),"Color",colour );
                        legend_Arr{plot_counter}= sprintf("type 01 Order %d",Order );
                        plot_counter = plot_counter+1;
                        hold on
                   case 5 
                      colour = colors{2};
                        colour(4) = alphas(j);
                        loglog( simul.beta_arr, abs(cell2mat(plot_structure(5,:))),"LineStyle", line_spec(j),"Color",colour    );
                        legend_Arr{plot_counter}= sprintf("type 05 Order %d",Order );
                        plot_counter = plot_counter+1;
                        hold on
                   otherwise
                       error("unknown type")
               end
        end  
         
        title( simul.title   )
        xlabel('$  \beta \cdot J$','Interpreter','latex')
        ylabel("relative err")
        legend(legend_Arr,'Location','northwest')
  
        figure(gcf)
        
    end  


    hold off
             
end

