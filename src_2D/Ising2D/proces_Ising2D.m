%notes: http://alps.comp-phys.org/mediawiki/index.php/Main_Page
close all
%sizes: 20*1.2.^(1:9)




ln_lopts = struct('Display', 0, 'maxit', 1);

filterPoints = 1;
%[names,X_crit_arr] = ising_names(3);
%skip = [1, 1, 0, 1];


[names,X_crit_arr] = ising_names(2);
skip = [1,1,3];

tbound = 0.08;

close all


critopts.T.name = "\beta";
critopts.T.crit = 1/8;

critopts.g.name = "\beta";
critopts.g.crit = 1/8;

plot_opts.marker_size = 2;

for j = 1:numel(names)

   
  
    
    
    X_crit = X_crit_arr(j);

    if skip(j) ~= 1
        
        fig = figure(j);
        set(fig, 'defaultAxesFontSize',16)
        
        t = tiledlayout(2,2 );
        
        
        %title(t,sprintf("2D transverse Ising, %s =%.1f", data.fixed_var, data.fixed_val));
        

        for i = 1:numel(names{j})

            data = fetch_matfiles(names{j}{i}, struct);
            data = filter_ising_results(data, struct('tol', 1e-10,'Tbound',[X_crit-  tbound,X_crit+  tbound ]));

            %d1=real(data.eps_i(:, 2)-data.eps_i(:, 1));
            %d2=real(  data.eps_i(:, 7) + data.eps_i(:, 6) + data.eps_i(:, 5) + data.eps_i(:, 4) + data.eps_i(:, 3)+ data.eps_i(:, 2) - 6 * data.eps_i(:, 1))/6;
            %d3=real(data.eps_i(:, 4) -  data.eps_i(:, 3));

            %plot(d1,d2,'*');

            %             marek_arr = 1/3 * real(  data.eps_i(:, 4) - data.eps_i(:, 1) +...
            %                 data.eps_i(:, 5) - data.eps_i(:, 2)+...
            %                 data.eps_i(:, 6)- data.eps_i(:, 3) );

            %             marek_arr = real((data.eps_i(:, 6) - data.eps_i(:, 1)) / 3 + ...
            %                 (data.eps_i(:, 5) - data.eps_i(:, 2)) / 2 + ...
            %                 (data.eps_i(:, 4) - data.eps_i(:, 3)));

            
            copts = critopts.(data.free_var);
            
            marek_arr = real(data.eps_i(:, 2) - data.eps_i(:, 1));

            plot_opts.colour = marek_arr;
            chi = data.chi(1);

            %T_c = 1.2737;  %https://journals.aps.org/prb/pdf/10.1103/PhysRevB.93.155157

            
             plot_m_vs_t(data, chi, marek_arr, i, j, X_crit, plot_opts,copts)
            
            plot_m_vs_t_marek(data, chi, marek_arr, i, j, X_crit, plot_opts,copts)
            
            plot_xi_marek(data, chi, marek_arr, i, j, X_crit, plot_opts,copts)
            
            plot_S_marek(data, chi, marek_arr, i, j, X_crit, plot_opts,copts)

        end
        
        x_width = 25;
        y_width = 20;
        
        set(fig, 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, x_width, y_width], 'PaperSize', [x_width, y_width])
      
        saveas(fig,'a.pdf')
    end
    
   
end

function name = plot_m_vs_t(data, chi, marek_arr, i, j, X_crit, plot_opts,copts)

    %figure(j);
    %subplot(2,2,1)

    nexttile(1)
    
    y_arr = data.m;
    x_arr = data.(data.free_var);

    if i == 1
        xline(X_crit, 'DisplayName', "$ T_c $")

        %title(sprintf("2D transverse Ising, %s =%.1f", data.fixed_var, data.fixed_val));

        switch data.free_var
            case 'T'
                xlabel("$T$", "Interpreter", "Latex");
            case 'g'
                xlabel("g", "Interpreter", "Latex");
        end

        ylabel("$\left < m \right >$", "Interpreter", "Latex");

        %legend('Location', 'westoutside', "Interpreter", "Latex")
        %lgd.Tile = 'east';
        %colorbar
        %set(gca, 'colorscale', 'log')
        %colormap('gray')
    end

    hold on

    plot(x_arr, y_arr, '*', 'MarkerSize', plot_opts.marker_size, 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    %scatter(x_arr, y_arr, 10, plot_opts.colour, 'filled', 'DisplayName', sprintf("$ \\chi  = %d $", chi))

    hold off

    name = sprintf("$ \\chi  = %d $", chi);
    
end

function plot_xi_marek(data, chi, marek_arr, i, j, X_crit, plot_opts,copts)

    %subplot(2,2,2)
    %nexttile
    nexttile(2)
    nu = 1;
    omega = 1;
    c = -0.5;
    d = 0.1;
    phi = 1;

    %y_arr =  (1./data.inv_corr_length .* marek_arr.^(1/nu))./( 1+c*marek_arr.^(omega)  )   ;
    %x_arr = (data.T - X_crit) .* (marek_arr.^(-1/1)) + d*marek_arr.^phi/nu;

    y_arr = log((1 ./ real(data.inv_corr_length) .* marek_arr.^(1 / nu)));
    x_arr = (data.(data.free_var) - X_crit) .* (marek_arr.^(-1/1));

    %figure(j + 600);

    if i == 1

        %title(sprintf("2D transverse Ising, %s =%.1f", data.fixed_var, data.fixed_val));

        switch data.free_var
            case 'T'
                xlabel("$ (T-T_c) \delta^{-1/\nu}$ ", "Interpreter", "Latex");
            case 'g'
                xlabel("g", "Interpreter", "Latex");
        end

        ylabel("$  -\log(\xi) - log(\delta)/ \nu  $", "Interpreter", "Latex");
        
       legend('Location', 'eastoutside', "Interpreter", "Latex") 
        
        %legend('Location', 'northwest', "Interpreter", "Latex")
        %colorbar
        %set(gca, 'colorscale', 'log')
        %colormap('gray')
    end

    hold on

    plot(x_arr, y_arr, '*', 'MarkerSize', plot_opts.marker_size, 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    %scatter(x_arr, y_arr, 10, plot_opts.colour, 'filled', 'DisplayName', sprintf("$ \\chi  = %d $", chi))

    hold off

end

function plot_S_marek(data, chi, marek_arr, i, j, X_crit, plot_opts,copts)
    %subplot(2,2,3)
    %nexttile
    c = 1/2;
nexttile(3)
    nu = 1;
    omega = 1.5;
    e = 0;
    d = 0;
    phi = 1;

    y_arr = log(exp(6 * data.S / c) .* marek_arr.^(1 / nu));
    x_arr = (data.(data.free_var) - X_crit) .* (marek_arr.^(-1/1)) + d * marek_arr.^phi / nu;

    %y_arr = 6*data.S/c +log(marek_arr) ;
    %x_arr =(data.T- X_crit).* (marek_arr.^(-1/1));

    %figure(j + 900);

    if i == 1
        %title(sprintf("2D transverse Ising, %s =%.1f", data.fixed_var, data.fixed_val));

        switch data.free_var
            case 'T'
                xlabel("$ (T-T_c) \delta^{-1/\nu} $", "Interpreter", "Latex");
            case 'g'
                xlabel("g", "Interpreter", "Latex");
        end
        
        ylabel("$  -\log(6S/c) - log(\delta)/ \nu  $", "Interpreter", "Latex");

        %legend('Location', 'northwest', "Interpreter", "Latex")
        %colorbar
        %set(gca, 'colorscale', 'log')
        %colormap('gray')
    end

    hold on

    plot(x_arr, y_arr, '*', 'MarkerSize', plot_opts.marker_size, 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    %scatter(x_arr, y_arr, 10, plot_opts.colour, 'filled', 'DisplayName', sprintf("$ \\chi  = %d $", chi))

end

function plot_m_vs_t_marek(data, chi, marek_arr, i, j, X_crit, plot_opts,copts)
    %subplot(2,2,4)
    nexttile(4)
    %nexttile
    nu = 1;
    c = 0;
    omega = 0;
    phi = 0;
    d = 0;

    yarr = data.m .* (marek_arr.^(- copts.crit )) ./ (1);
    %xarr = (data.T - X_crit) / X_crit .* (marek_arr.^(-1 / nu)) + d * marek_arr.^(phi / nu);

    xarr = (data.(data.free_var) - X_crit) / X_crit .* (marek_arr.^(-1 / nu)) + d * marek_arr.^(phi / nu);

    
    %figure(j + 400);

    if i == 1
        %title(sprintf("2D transverse Ising, %s =%.1f", data.fixed_var, data.fixed_val));

        switch data.free_var
            case 'T'
                xlabel("$ (T-T_c) \delta^{-1/\nu} $", "Interpreter", "Latex");
            case 'g'
                xlabel("g", "Interpreter", "Latex");
        end
        
        ylabel(  sprintf("$  m \\delta ^{-%s/\\nu}$  ",copts.name), "Interpreter", "Latex");

        %legend('Location', 'southwest', "Interpreter", "Latex")

        %set(gca, 'ColorScale', 'log')

        %colorbar
        %set(gca, 'colorscale', 'log')
        %colormap('gray')
    end

    hold on

    plot(xarr, yarr, '*', 'MarkerSize', plot_opts.marker_size, 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    %scatter(xarr, yarr, 10, plot_opts.colour, 'filled', 'DisplayName', sprintf("$ \\chi  = %d $", chi))

    hold off

end
