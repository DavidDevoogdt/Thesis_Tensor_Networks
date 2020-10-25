function test
%mpo_type_comparison_exact_generic_order;
%mpo_type_comparison_M
phase_transition();
end



function phase_transition()


%change this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_z_2 = [1,0;0,-1];
S_x_2 = [0,1;1,0];

generate_opts.testing=0;
opt3 = struct([]);
opt5.method="diag";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta_max = 20; %from simulations
t_max = 3;

simul_struct.order = 3;
simul_struct.mpo_type = 4;
simul_struct.observable = 'A';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch simul_struct.observable
    case 'sx'
        simul_struct.title = 'S_x';
        simul_struct.Z = 'S_x';
    case 'xi'
        simul_struct.title = 'xi';
        simul_struct.Z = 'xi';
    case 'A'
        simul_struct.title = 'A';
        simul_struct.Z = 'A';
    otherwise
        error('unknown observable')
end


t_arr = exp( log(1/beta_max):0.2:log(t_max));
beta_arr = 1./t_arr;

%beta_arr = exp( l:0.1:log(50));
%t_arr = 1./beta_arr;

beta_len = size(beta_arr,2);

gamma_arr = 0.2:0.1:3;
gamma_len = size(gamma_arr,2);

results = zeros(beta_len,gamma_len);

disp(gamma_arr);

for beta_i = 1:beta_len
    beta = beta_arr(beta_i);
    fprintf("beta: %.4f gamma: ",beta)
    for gamma_i = 1:gamma_len
        gamma = gamma_arr(gamma_i);
        model_opts.g = gamma;
        [simul,H_1_tensor,H_2_tensor,opt4,d] = models('t_ising',model_opts);
        
        
        switch simul_struct.mpo_type
            case 3           
                MPO = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,3,simul_struct.order,opt3,generate_opts);
            case 4
                MPO = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,4,simul_struct.order,opt4,generate_opts);
            otherwise
                error('unknown mpo type')
        end
        
        
        MPO_matrix = MPO.MPO_cell * MPO.nf;
        
        switch simul_struct.observable
            case 'sx'
                m = get_expectation(MPO_matrix, S_x_2);
                results(beta_i,gamma_i) =  m ;
            case 'xi'
                [A, xi_inv] =   get_correlation_length(MPO_matrix, S_z_2,2);
                fprintf(" %.4e ",xi_inv);
                results(beta_i,gamma_i) = 1/ xi_inv ;
            case 'A'
                [A, xi_inv] =   get_correlation_length(MPO_matrix, S_z_2,2);
                fprintf(" %.4e ",A);
                results(beta_i,gamma_i) = A ;
            otherwise
                error('unknown observable')
        end
        
    end
    fprintf("\n")
    
end

figure()

h=gca;

s = surf(gamma_arr-1,t_arr,results);
colorbar
xlabel('h')
ylabel('T')
zlabel(simul_struct.Z);
title(simul_struct.title);

s.EdgeColor = 'none';
if (simul_struct.observable == 'xi') || (simul_struct.observable == 'A')
    set(h,'zscale','log')
    set(gca,'ColorScale','log')
      
end
end




function mpo_type_comparison_exact_generic_order

%change this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

generate_opts.testing=0;

opt3 = struct([]);


model_opts.g =  0.5;
[simul,H_1_tensor,H_2_tensor,opt4,d] = models('t_ising',model_opts);


simul.Order_arr = [3,4,5,6];
simul.types = [3,4];
simul.M = 9;
simul.beta_arr = 10.^(  -3:0.1:log10(50) );
simul.cyclic = 0;


opt5.method="diag";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%plot stuff
order_size = size(simul.Order_arr,2);
line_spec = ["-"    ,"--"   ,"-."   ,":"    ,"-"    ,"--"   ,"-."   ,":"];
alphas= [1     ,1      ,1      ,1      ,0.5    ,0.5    ,0.5    ,0.5];
colors = {[1 0 0 0],[0 1 0 0],[0 0 1 0]}; %reserved for specific type

legend_Arr = cell(size(simul.types,2)*order_size,1);
legend_Arr(:) = {"todo"};
beta_len = size(simul.beta_arr,2);

plot_counter = 1;
%hold off
figure();
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
compare_opts.cyclic=simul.cyclic;

%loop over different orders
for j = 1:order_size
    Order = simul.Order_arr(j);
    plot_structure = cell( 5,beta_len);
    
    
    
    %loop over temps
    for i = 1:beta_len
        beta = simul.beta_arr(i);
        mpo_base = MPO_beta_N{i};
        mpo_base_matrix = mpo_base.H_exp(simul.M-1,1,simul.cyclic); %buffered, not calculted again every time
        
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
    legend(legend_Arr,'Location','northwest','NumColumns',2)
    
    figure(gcf)
    
end


hold off

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mpo_type_comparison_M

%change this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

generate_opts.testing=0;

opt3 = struct([]);

[simul,H_1_tensor,H_2_tensor,opt4] = models('t_ising');
simul.Order_arr = [2,3,4,5,6];
simul.types = [3,4];
simul.M = 10;
simul.beta_arr = 10.^(  -3:0.1:log10(50) );
simul.cyclic = 1;


opt5.method="diag";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%plot stuff
M_size = size(simul.M_arr,2);
line_spec = ["-"    ,"--"   ,"-."   ,":"    ,"-"    ,"--"   ,"-."   ,":"];
alphas= [1     ,1      ,1      ,1      ,0.5    ,0.5    ,0.5    ,0.5];
colors = {[1 0 0 0],[0 1 0 0],[0 0 1 0]}; %reserved for specific type

legend_Arr = cell(size(simul.types,2)*M_size,1);
legend_Arr(:) = {"todo"};
beta_len = size(simul.beta_arr,2);

plot_counter = 1;
%hold off
figure();
%%%


%
compare_opts.ref=2;
compare_opts.cyclic=simul.cyclic;


plot_structure = cell( 5,beta_len,M_size);


%loop over temps
for i = 1:beta_len
    beta = simul.beta_arr(i);
    
    mpo_base =   generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,0);
    
    for t=1:size(simul.types,2)
        switch simul.types(t)
            case 3
                mpo_3 = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,3,simul.Order,opt3,generate_opts);
            case 4
                mpo_4 = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,4,simul.Order,opt4,generate_opts);
            case 5
                mpo_5 = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,5,simul.Order,opt5,generate_opts);
            otherwise
                error("unknown type")
        end
    end
    


    for j = 1:M_size
        M = simul.M_arr(j);
        
        
        mpo_base_matrix = mpo_base.H_exp(M-1,1,simul.cyclic);

        fprintf("beta %.4e M %d",beta,M);
        for t=1:size(simul.types,2)
            switch simul.types(t)
                case 3
                    err_03 = error_eigenvalue(mpo_3,mpo_base_matrix,"array",M,d,compare_opts);
                    fprintf(" err 03 %.4e",err_03);
                    plot_structure{3,i,j} = err_03;

                case 4
                    err_04 = error_eigenvalue(mpo_4,mpo_base_matrix,"array",M,d,compare_opts);
                    fprintf(" err 04 %.4e",err_04);
                    plot_structure{4,i,j} = err_04;
                case 5
                    err_05 = error_eigenvalue(mpo_5,mpo_base_matrix,"array",M,d,compare_opts);
                    fprintf(" err 05 %.4e",err_05);
                    plot_structure{5,i,j} = err_05;
                otherwise
                    error("unknown type")
            end
        end
        fprintf("\n");

    end

end

%plotting loop for current order
for j=1:M_size
    M = simul.M_arr(j);
    for t=1:size(simul.types,2)
        switch simul.types(t)
            case 3
                colour = colors{1};
                colour(4) = alphas(j);
                loglog( simul.beta_arr, abs(cell2mat(plot_structure(3,:,j))), "LineStyle", line_spec(j),"Color",colour );
                legend_Arr{plot_counter}= sprintf("type 03 M %d",M );
                plot_counter = plot_counter+1;
                hold on
            case 4
                colour = colors{3};
                colour(4) = alphas(j);
                loglog( simul.beta_arr, abs(cell2mat(plot_structure(4,:,j))),"LineStyle", line_spec(j),"Color",colour );
                legend_Arr{plot_counter}= sprintf("type 01  M %d",M );
                plot_counter = plot_counter+1;
                hold on
            case 5
                colour = colors{2};
                colour(4) = alphas(j);
                loglog( simul.beta_arr, abs(cell2mat(plot_structure(5,:,j))),"LineStyle", line_spec(j),"Color",colour    );
                legend_Arr{plot_counter}= sprintf("type 05  M %d",M );
                plot_counter = plot_counter+1;
                hold on
            otherwise
                error("unknown type")
        end
    end
end
title( simul.title   )
xlabel('$  \beta \cdot J$','Interpreter','latex')
ylabel("relative error")
legend(legend_Arr,'Location','northwest','NumColumns',2)

figure(gcf)

hold off

end


