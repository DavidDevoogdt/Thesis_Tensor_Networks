function test
mpo_type_comparison_exact_generic_order;
%mpo_type_comparison_M
%phase_transition();

end



function phase_transition()


%change this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_z_2 = [1,0;0,-1];
S_x_2 = [0,1;1,0];

generate_opts.testing=0;
generate_opts.MPO_type="matrix";

opt3 = struct([]);
opt5.method="svd";

a=0.5;
model_opts.J = 1/(2*a);

square_size = 4;
step_size = 0.3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %from simulations


simul_struct.order = 5;
simul_struct.mpo_type = 4;
simul_struct.observable = 'xi';


simul_struct.cor_min = 1;
simul_struct.cor_max = 4;

simul_struct.model = 't_ising';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta_max = 10; %saveguard for numerical errors

switch simul_struct.observable
    case 'sx'
        simul_struct.title = 'S_x';
        simul_struct.Z = 'S_x';
        
        t_max = 8;

        t_arr = exp( log(1/beta_max):0.2:log(t_max));
        %t_arr =  1/beta_max :0.5:t_max;   %exp( log(1/beta_max):0.2:log(t_max))
        beta_arr = 1./t_arr;
        
        %beta_arr = exp( l:0.1:log(50));
        %t_arr = 1./beta_arr;
        
        beta_len = size(beta_arr,2);
        
        gamma_arr = 0.2:0.1:3;
        gamma_len = size(gamma_arr,2);
        
        
    case 'xi'
        simul_struct.title = '1/xi';
        simul_struct.title2 = 'A';
        simul_struct.Z = '1/xi';
        simul_struct.Z2 = 'A';

        m_arr = -square_size:step_size:square_size;
        T_arr_max = square_size;

        gamma_arr = 1-m_arr*a
        gamma_arr = gamma_arr( gamma_arr>0.4);

        gamma_len = size(gamma_arr,2);


        t_min =  model_opts.J/beta_max;

        t_arr =  t_min:step_size:T_arr_max
        beta_arr = 1./t_arr;
        beta_len = size(beta_arr,2);
        

    otherwise
        error('unknown observable')
end




% 
 results = zeros(beta_len,gamma_len,2);


disp(gamma_arr);

for beta_i = 1:beta_len
    beta = beta_arr(beta_i);
    fprintf("beta: %.4f gamma: ",beta)
    for gamma_i = 1:gamma_len
        gamma = gamma_arr(gamma_i);
        
        
        model_opts.g = gamma;
        
        [simul,H_1_tensor,H_2_tensor,opt4,d] = models(simul_struct.model,model_opts);
        
        
        switch simul_struct.mpo_type
            case 2
                MPO = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,2,simul_struct.order,opt3,generate_opts);
            case 3           
                MPO = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,3,simul_struct.order,opt3,generate_opts);
            case 4
                MPO = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,4,simul_struct.order,opt4,generate_opts);
            otherwise
                error('unknown mpo type')
        end
        
        %MPO.nf
        
        MPO_matrix = MPO.MPO_cell;
        
        switch simul_struct.observable
            case 'sx'
                m = get_expectation(MPO_matrix, S_x_2,MPO.nf);
                 fprintf(" %.4e ",m);
                 fprintf(" %f ",gamma);
                results(beta_i,gamma_i,1) =  m ;
            case 'xi'
                [A, xi_inv] =   get_correlation_length(MPO_matrix, S_z_2,simul_struct.cor_min,simul_struct.cor_max,a,MPO.nf);
                fprintf(" %.4e ",xi_inv);
                results(beta_i,gamma_i,1) = xi_inv ;
                results(beta_i,gamma_i,2) = A ;
            otherwise
                error('unknown observable')
        end
        
    end
    fprintf("\n")
    
end

%%

% a = 0.1;
% J= 1/(2*a);

switch simul_struct.observable 
    case 'xi'

    Z= model_opts.J^(-1/4);

    T = t_arr;

    m= (1-gamma_arr)/a;

    [m_grid,T_grid] = meshgrid(m,T);

    region_I = (T_grid < m_grid) & (m_grid>0);
    region_II = (T_grid < -m_grid) & (m_grid<0);
    region_III = ~(region_I | region_II );

    m_I = m_grid; T_I = T_grid;
    m_I(~region_I) = nan; 
    T_I(~region_I) = nan; 

    m_II = m_grid; T_II = T_grid;
    m_II(~region_II) = nan; 
    T_II(~region_II) = nan; 

    m_III = m_grid; T_III = T_grid;
    m_III(~region_III) = nan; 
    T_III(~region_III) = nan; 



    xi_inv_I = (  (2*m_I.*T_I/3.1415).^(0.5).*exp(-m_I./T_I));
    xi_inv_III = (3.1415*T_III/4);
    xi_inv_II = abs(m_II) + (2*abs(m_II).*T_II/3.1415).^(0.5).*exp(-abs(m_II)./T_II);


    figure();
    h=gca;

    s = surf(m_grid,T_grid,results(:,:,1));

    colorbar
    xlabel('m')
    ylabel('T')
    zlabel(simul_struct.Z);
    title(simul_struct.title);
    s.EdgeColor = 'none';

    hold on;
    surf(m_I,T_I,xi_inv_I,'FaceAlpha',0.5);
    surf(m_II,T_II,xi_inv_II,'FaceAlpha',0.5);
    surf(m_III,T_III,xi_inv_III,'FaceAlpha',0.5);



    figure();
    h=gca;



    A_prefact = Z* (T_grid).^(0.25);

    A_I = A_prefact.*(m_I./T_I).^(0.25);
    A_III = A_prefact.*0.85871456;
    A_II = A_prefact.*abs(m_II./T_II).^(-0.75);

    s = surf(m_grid,T_grid,results(:,:,2));

    colorbar
    xlabel('m')
    ylabel('T')
    zlabel(simul_struct.Z2);
    title(simul_struct.title2);

    s.EdgeColor = 'none';

    hold on;
    surf(m_I,T_I,A_I,'FaceAlpha',0.5);
    surf(m_II,T_II, A_II,'FaceAlpha',0.5);
    surf(m_III,T_III,A_III,'FaceAlpha',0.5);

    hold off;
    case 'sx'
        figure();
        

        s = surf(gamma_arr,t_arr,real(results(:,:,1)));

        colorbar
        xlabel('gamma')
        ylabel('T')
        zlabel(simul_struct.Z);
        title(simul_struct.title);
        s.EdgeColor = 'none';
    otherwise
        error('unknown plot ype')
end

end

%surf( m,t_arr,xi_inv_calc );


%%
function mpo_type_comparison_exact_generic_order

simulatiemodellen = ["t_ising","random","random","random","Heisenberg_2D","Heisenberg_2D_X","random","random"];
%simulatiemodellen = ["Heisenberg_2D"];
models_len = size(simulatiemodellen,2);

    for round = 1:models_len
%change this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    generate_opts.testing=0;
    generate_opts.MPO_type = "matrix";

    opt3 = struct([]);


    model_opts.g =  0.5;
    %[simul,H_1_tensor,H_2_tensor,opt4,d] = models('Heisenberg_2D',model_opts);
    %[simul,H_1_tensor,H_2_tensor,opt4,d] = models('t_ising',model_opts);
    [simul,H_1_tensor,H_2_tensor,opt4,d] = models( simulatiemodellen(round)  ,model_opts);

    opt4.single_threshold = 1e-12;
%     
% 
% %simulatiemodellen = ["random","random","random","Heisenberg_2D","t_ising","Heisenberg_2D_X","random","random"];
% simulatiemodellen = [1e-10,1e-11,1e-12,1e-13,1e-14];
% models_len = size(simulatiemodellen,2);
% 
% generate_opts.testing=0;
% generate_opts.MPO_type = "matrix";
% 
% opt3 = struct([]);
% 
% 
% model_opts.g =  0.5;
% %[simul,H_1_tensor,H_2_tensor,opt4,d] = models('Heisenberg_2D',model_opts);
% %[simul,H_1_tensor,H_2_tensor,opt4,d] = models('t_ising',model_opts);
% [simul,H_1_tensor,H_2_tensor,opt4,d] = models( "t_ising"  ,model_opts);
% 
%     for round = 1:models_len
% 
%     opt4.single_threshold =simulatiemodellen(round);
%     


    simul.Order_arr = [2,4];
    simul.types = [2,5];
    simul.M =8 ;
    simul.beta_arr = 10.^(  -3:0.1:1.5);
    %simul.beta_arr = 10.^(  -0:0.05:1);
    simul.cyclic = 1;


    opt5.method="svd";



    compare_opts.ref=2;
    compare_opts.cyclic=simul.cyclic;

    opt2 = {};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%plot stuff
    order_size = size(simul.Order_arr,2);
    line_spec = ["-"    ,"--"   ,"-."   ,":"    ,"-"    ,"--"   ,"-."   ,":"];
    alphas= [1     ,1      ,1      ,1      ,0.5    ,0.5    ,0.5    ,0.5];
    colors = {[1 0 0 0],[0 1 0 0],[0 0 1 0],[0,1,1]}; %reserved for specific type

    legend_Arr = cell(size(simul.types,2)*order_size,1);
    legend_Arr(:) = {"todo"};
    beta_len = size(simul.beta_arr,2);

    plot_counter = 1;
    %hold off
    figure();
    x_width = 15;
    y_width=10;
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0,0,x_width,y_width],'PaperSize', [x_width,y_width])
    %%%


    %pregenerate all mpo's
    MPO_beta_N = cell(beta_len,1);
    for i=1:beta_len
        beta = simul.beta_arr(i);

        %generateMPO(d,H_1_tensor,H_2_tensor,type,order,type_opts,opts)
        MPO_beta_N{i} = generateMPO(simul.d,-beta*H_1_tensor,-beta*H_2_tensor,0 );
    end

    %


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
                    case 2   
                        mpo_2 = generateMPO(d,-beta*H_1_tensor,-beta*H_2_tensor,2,Order,opt2,generate_opts);
                        err_02 = error_eigenvalue(mpo_2,mpo_base_matrix,"array",simul.M,d,compare_opts);
                        fprintf(" err 02 %.4e",err_02);
                        plot_structure{2,i} = err_02;
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
                case 2
                    colour = colors{1};
                    colour(4) = alphas(j);
                    loglog( simul.beta_arr, abs(cell2mat(plot_structure(2,:))), "LineStyle", line_spec(j),"Color",colour );
                    legend_Arr{plot_counter}= sprintf("B:%d",Order );
                    plot_counter = plot_counter+1;
                    hold on
                 case 4
                    colour = colors{3};
                    colour(4) = alphas(j);
                    loglog( simul.beta_arr, abs(cell2mat(plot_structure(4,:))),"LineStyle", line_spec(j),"Color",colour );
                    legend_Arr{plot_counter}= sprintf("A:%d",Order );
                    plot_counter = plot_counter+1;
                    hold on    
                case 3
                    colour = colors{2};
                    colour(4) = alphas(j);
                    loglog( simul.beta_arr, abs(cell2mat(plot_structure(3,:))), "LineStyle", line_spec(j),"Color",colour );
                    legend_Arr{plot_counter}= sprintf("C:%d",Order );
                    plot_counter = plot_counter+1;
                    hold on
                case 5
                    colour = colors{4};
                    colour(4) = alphas(j);
                    loglog( simul.beta_arr, abs(cell2mat(plot_structure(5,:))),"LineStyle", line_spec(j),"Color",colour    );
                    legend_Arr{plot_counter}= sprintf("D:%d",Order );
                    plot_counter = plot_counter+1;
                    hold on
                otherwise
                    error("unknown type")
            end
        end

        title( simul.title   )
        xlabel('$  \beta \cdot J$','Interpreter','latex', 'FontSize', 11)
        ylabel('$  \epsilon $','Interpreter','latex', 'FontSize', 11)
        legend(legend_Arr,'Location','northwest','NumColumns',2 )

        ylim([0 10])

        figure(gcf)

    end


    hold off

    filename = sprintf('../auto_fig/comp%s.pdf',datestr(now,'mm-dd-yy_HH-MM-SS'));
    saveas(gcf,filename )

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mpo_type_comparison_M

%change this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

generate_opts.testing=0;
generate_opts.MPO_type="cell";

opt3 = struct([]);

model_opts.g =  0.5;

[simul,H_1_tensor,H_2_tensor,opt4] = models('Heisenberg_2D',model_opts);
simul.M_arr = [6,7,8,9];
simul.types = [3,4];
simul.Order = 4;
simul.beta_arr = 10.^(  -3.1:0.1:log10(20) );
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
x_width = 15;
y_width=10;
set(gcf,'PaperUnits','centimeters','PaperPosition',[0,0,x_width,y_width],'PaperSize', [x_width,y_width])
%%%


%
compare_opts.ref=2;
compare_opts.cyclic=simul.cyclic;


plot_structure = cell( 5,beta_len,M_size);


%loop over temps
for i = 1:beta_len
    beta = simul.beta_arr(i);
    
    mpo_base =   generateMPO(simul.d,-beta*H_1_tensor,-beta*H_2_tensor,0);
    
    for t=1:size(simul.types,2)
        switch simul.types(t)
            case 3
                mpo_3 = generateMPO(simul.d,-beta*H_1_tensor,-beta*H_2_tensor,3,simul.Order,opt3,generate_opts);
            case 4
                mpo_4 = generateMPO(simul.d,-beta*H_1_tensor,-beta*H_2_tensor,4,simul.Order,opt4,generate_opts);
            case 5
                mpo_5 = generateMPO(simul.d,-beta*H_1_tensor,-beta*H_2_tensor,5,simul.Order,opt5,generate_opts);
            otherwise
                error("unknown type")
        end
    end
    
    fprintf("beta %.4e \n",beta);

    for j = 1:M_size
        M = simul.M_arr(j);

        mpo_base_matrix = mpo_base.H_exp(M-1,1,simul.cyclic);

        fprintf(" M %d ",M);
        for t=1:size(simul.types,2)
            switch simul.types(t)
                case 3
                    err_03 = error_eigenvalue(mpo_3,mpo_base_matrix,"array",M,simul.d,compare_opts);
                    fprintf(" err 03 %.4e",err_03);
                    plot_structure{3,i,j} = err_03;

                case 4
                    err_04 = error_eigenvalue(mpo_4,mpo_base_matrix,"array",M,simul.d,compare_opts);
                    fprintf(" err 04 %.4e",err_04);
                    plot_structure{4,i,j} = err_04;
                case 5
                    err_05 = error_eigenvalue(mpo_5,mpo_base_matrix,"array",M,simul.d,compare_opts);
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
                legend_Arr{plot_counter}= sprintf("B-%d",M );
                plot_counter = plot_counter+1;
                hold on
            case 4
                colour = colors{3};
                colour(4) = alphas(j);
                loglog( simul.beta_arr, abs(cell2mat(plot_structure(4,:,j))),"LineStyle", line_spec(j),"Color",colour );
                legend_Arr{plot_counter}= sprintf("A-%d",M );
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
ylabel('$ \epsilon $','Interpreter','latex')
legend(legend_Arr,'Location','northwest','NumColumns',2)

figure(gcf)

hold off


filename = sprintf('../auto_fig/compM%s.pdf',datestr(now,'mm-dd-yy_HH-MM-SS'));
saveas(gcf,filename )



end


