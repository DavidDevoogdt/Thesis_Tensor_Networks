%notes: http://alps.comp-phys.org/mediawiki/index.php/Main_Page

%sizes: 20*1.2.^(1:9)

names = {{
%     'Ising2D_g=2.5000e+00_chi=20_22_March_2021_10:11_a';
%     'Ising2D_g=2.5000e+00_chi=25_22_March_2021_10:11_a';
%     'Ising2D_g=2.5000e+00_chi=30_22_March_2021_10:11_a';
%     'Ising2D_g=2.5000e+00_chi=40_22_March_2021_10:46_a';
%     'Ising2D_g=2.5000e+00_chi=60_22_March_2021_11:49_a';
%     'Ising2D_g=2.5000e+00_chi=65_22_March_2021_13:52_a';
%     'Ising2D_g=2.5000e+00_chi=70_22_March_2021_13:51_a';
%     },{
%     'Ising2D_g=2.5000e+00_chi=8_02_April_2021_09:38';
%     'Ising2D_g=2.5000e+00_chi=11_02_April_2021_09:38';
%     'Ising2D_g=2.5000e+00_chi=16_02_April_2021_09:38';
%     'Ising2D_g=2.5000e+00_chi=23_02_April_2021_09:38';
%     'Ising2D_g=2.5000e+00_chi=32_02_April_2021_09:38';
%     'Ising2D_g=2.5000e+00_chi=45_02_April_2021_09:38';
%     },{
%     'Ising2D_g=1.5000e+00_chi=8_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=11_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=16_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=23_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=32_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=45_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=64_02_April_2021_11:17';
%     },{
%     'Ising2D_g=0.0000e+00_chi=8_02_April_2021_11:46';
%     'Ising2D_g=0.0000e+00_chi=11_02_April_2021_11:46';
%     'Ising2D_g=0.0000e+00_chi=16_02_April_2021_11:46';
%     'Ising2D_g=0.0000e+00_chi=23_02_April_2021_11:46';
%     'Ising2D_g=0.0000e+00_chi=32_02_April_2021_11:46';
%     'Ising2D_g=0.0000e+00_chi=45_02_April_2021_11:46';
%     },{
    'Ising2D_g=0.0000e+00_chi=8_08_April_2021_09:52';
    'Ising2D_g=0.0000e+00_chi=11_08_April_2021_10:11';
    'Ising2D_g=0.0000e+00_chi=16_08_April_2021_10:43';
    'Ising2D_g=0.0000e+00_chi=23_08_April_2021_11:50';
    'Ising2D_g=0.0000e+00_chi=32_08_April_2021_14:41';
    'Ising2D_g=0.0000e+00_chi=45_08_April_2021_18:59';
    'Ising2D_g=0.0000e+00_chi=64_09_April_2021_05:59';
    },{
    'Ising2D_g=1.5000e+00_chi=8_08_April_2021_10:15';
    'Ising2D_g=1.5000e+00_chi=11_08_April_2021_10:28';
    'Ising2D_g=1.5000e+00_chi=16_08_April_2021_10:48';
    'Ising2D_g=1.5000e+00_chi=23_08_April_2021_11:29';
    'Ising2D_g=1.5000e+00_chi=32_08_April_2021_12:19';
    'Ising2D_g=1.5000e+00_chi=45_08_April_2021_16:50';
    'Ising2D_g=1.5000e+00_chi=64_09_April_2021_04:35';
    },{
    'Ising2D_g=2.5000e+00_chi=8_08_April_2021_09:52';
    'Ising2D_g=2.5000e+00_chi=11_08_April_2021_10:02';
    'Ising2D_g=2.5000e+00_chi=16_08_April_2021_10:17';
    'Ising2D_g=2.5000e+00_chi=23_08_April_2021_10:52';
    'Ising2D_g=2.5000e+00_chi=32_08_April_2021_12:39';
    'Ising2D_g=2.5000e+00_chi=45_08_April_2021_17:46';
    'Ising2D_g=2.5000e+00_chi=64_09_April_2021_05:52';
    },{
    'Ising2D_g=2.9000e+00_chi=8_08_April_2021_09:52';
    'Ising2D_g=2.9000e+00_chi=11_08_April_2021_10:03';
    'Ising2D_g=2.9000e+00_chi=16_08_April_2021_10:06';
}};


%1.27375%1.2737(2)

T_c_arr = [  2/log( 1+ sqrt(2) ),1.980,1.27348, 0.9  ];

ln_lopts = struct('Display', 0, 'maxit', 1);

filterPoints = 1;


skip = [1,1,0,1];

close all

plot_opts.marker_size = 2;


for j = 1:numel(names)
    
    T_c = T_c_arr(j);
    
    if skip(j)~=1

        for i = 1:numel(names{j})

            data = fetch_matfiles(names{j}{i}, struct);
            data = filter_ising_results(data, struct);

            %d1=real(data.eps_i(:, 2)-data.eps_i(:, 1));
            %d2=real(  data.eps_i(:, 7) + data.eps_i(:, 6) + data.eps_i(:, 5) + data.eps_i(:, 4) + data.eps_i(:, 3)+ data.eps_i(:, 2) - 6 * data.eps_i(:, 1))/6;
            %d3=real(data.eps_i(:, 4) -  data.eps_i(:, 3));
            
            
            %plot(d1,d2,'*');

            
%             marek_arr = 1/3 * real(  data.eps_i(:, 4) - data.eps_i(:, 1) +...
%                 data.eps_i(:, 5) - data.eps_i(:, 2)+...
%                 data.eps_i(:, 6)- data.eps_i(:, 3) );
            
           marek_arr = real((data.eps_i(:, 6) - data.eps_i(:, 1)) / 3 + ...
                        (data.eps_i(:, 5) - data.eps_i(:, 2)) / 2 + ...
                        (data.eps_i(:, 4) - data.eps_i(:, 3)));
%                       
%            marek_arr =real(  data.eps_i(:, 2) -  data.eps_i(:, 1));

            
            
            chi = data.chi(1);

            %T_c = 1.2737;  %https://journals.aps.org/prb/pdf/10.1103/PhysRevB.93.155157

            plot_m_vs_t(data, chi,marek_arr, i, j, T_c)
            plot_m_vs_t_marek(data, chi,marek_arr, i, j, T_c,plot_opts)


            plot_xi_marek(data, chi,marek_arr, i, j, T_c,plot_opts)

            plot_S_marek(data, chi,marek_arr, i, j, T_c,plot_opts)

        end

    end
end

function plot_m_vs_t(data, chi,marek_arr, i, j, T_c)

    figure(j);
    
    
    

    if i == 1
        xline( T_c , 'DisplayName', "$ T_c $")

        title(sprintf("2D transverse Ising, g=%.4f", data.model_params.g));
        xlabel("$\frac{k T}{J}$", "Interpreter", "Latex", 'FontSize', 12);
        ylabel("$\left < m \right >$", "Interpreter", "Latex");

        legend('Location', 'southwest', "Interpreter", "Latex", 'FontSize', 12)
    end
    
    hold on
    plot(data.T, data.m, '*-', 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    
    hold off

end


function plot_xi_marek(data, chi,marek_arr, i, j, T_c,plot_opts)

    nu=1;
    omega=1;
    c=-0.5;
    d=0.1;
    phi=1;
    


    %y_arr =  (1./data.inv_corr_length .* marek_arr.^(1/nu))./( 1+c*marek_arr.^(omega)  )   ; 
    %x_arr = (data.T - T_c) .* (marek_arr.^(-1/1)) + d*marek_arr.^phi/nu;

    y_arr =  (1./data.inv_corr_length .* marek_arr.^(1/nu))   ; 
    x_arr = (data.T - T_c) .* (marek_arr.^(-1/1));
    
    figure(j + 600);

    if i == 1

        title(sprintf("$$H =ZZ+ %.2f X$$, $$T_c = %.5f$$ ", data.model_params.g,T_c),'interpreter','latex');
        xlabel("$(T-T_c)\delta^{-1/\nu}$", "Interpreter", "Latex", 'FontSize', 12);
        ylabel("$ \xi \delta  $", "Interpreter", "Latex");

         legend('Location', 'northwest', "Interpreter", "Latex", 'FontSize', 12)

    end

    hold on
    plot(x_arr, real(y_arr), '*', 'MarkerSize',plot_opts.marker_size,'DisplayName', sprintf("$ \\chi  = %d $", chi));
    %scatter(x_arr, real(y_arr), 10, marek_arr, 'filled','DisplayName', sprintf("$ \\chi  = %d $", chi))

    hold off

end

function plot_S_marek(data, chi,marek_arr, i, j, T_c,plot_opts)
    c = 1/2;

    nu=1;
    omega=1.5;
    e=0;
    d=0;
    phi=1;
    


    y_arr = log(  exp(6*data.S/c) .* marek_arr.^(1/nu)./( 1+e*marek_arr.^omega )  )   ; 
    x_arr = (data.T - T_c) .* (marek_arr.^(-1/1)) + d*marek_arr.^phi/nu;
    
    
    %y_arr = 6*data.S/c +log(marek_arr) ;
    %x_arr =(data.T- T_c).* (marek_arr.^(-1/1));

    figure(j + 900);

    if i == 1

        title(sprintf("$$H =ZZ+ %.2f X$$, $$T_c = %.5f$$ ", data.model_params.g,T_c),'interpreter','latex');
        xlabel("$(T-T_c)  \delta^{-1/\nu}$", "Interpreter", "Latex", 'FontSize', 12);
        ylabel("$ 6 S / c + \ln \delta  $", "Interpreter", "Latex");

        legend('Location', 'northwest', "Interpreter", "Latex", 'FontSize', 12)

        
        %colorbar
    end

    hold on
    plot(x_arr, real(y_arr), '*', 'MarkerSize',plot_opts.marker_size,'DisplayName', sprintf("$ \\chi  = %d $", chi));
    %scatter(x_arr, real(y_arr), 10, marek_arr, 'filled','DisplayName', sprintf("$ \\chi  = %d $", chi))

    hold off

end

function plot_m_vs_t_marek(data, chi,marek_arr, i, j, T_c,plot_opts)

    nu=1;
    c=-0.3;
    omega = 1.4;
    phi=1;
    d=-0.1;
    
    yarr = data.m .* (marek_arr.^(-1/8)) ./( 1 + c* marek_arr.^(omega) );
    xarr = (data.T - T_c) / T_c .* (marek_arr.^(-1/nu))+d*marek_arr.^(phi/nu);

  
    figure(j + 400);

    if i == 1
        title(sprintf("$$H =ZZ+ %.2f X$$, $$T_c = %.5f$$ ", data.model_params.g,T_c),'interpreter','latex');
        xlabel("$ (T-T_C) \delta^{-1 / \nu} $", "Interpreter", "Latex", 'FontSize', 12);
        ylabel("$  \left< m \right>  \delta ^{ -\beta / \nu}  $", "Interpreter", "Latex");

        legend('Location', 'southwest', "Interpreter", "Latex", 'FontSize', 12)

        %set(gca, 'ColorScale', 'log')

        %colorbar
        %set(gca, 'colorscale', 'log')
        %colormap('gray')
    end

    hold on

    plot(xarr, yarr,  '*', 'MarkerSize',plot_opts.marker_size, 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    %scatter(xarr, yarr, 10, marek_arr, 'filled','DisplayName', sprintf("$ \\chi  = %d $", chi))
   

    hold off

end
