%notes: http://alps.comp-phys.org/mediawiki/index.php/Main_Page

%sizes: 20*1.2.^(1:9)

%[names,T_c_arr] = ising_names(2);

names = {{%'Ising2D_T=5.0000e-01_chi=8_26_April_2021_11:28';
%     'Ising2D_T=5.0000e-01_chi=12_26_April_2021_11:34';
%     },{
%     'Ising2D_g=2.5000e+00_chi=8_26_April_2021_14:30';
%     },{
%'Ising2D_g=2.5000e+00_chi=30_02_May_2021_10:27',
%'Ising2D_g=2.5000e+00_chi=9_09_May_2021_12:03';
%'Ising2D_g=2.5000e+00_chi=16_09_May_2021_12:03',
%'Ising2D_g=2.5000e+00_chi=8_09_May_2021_14:27',
%'Ising2D_g=2.5000e+00_chi=8_09_May_2021_15:30',
%'Ising2D_g=2.5000e+00_chi=8_09_May_2021_15:55',
%'Ising2D_g=2.5000e+00_chi=8_09_May_2021_16:01';
%'Ising2D_g=2.5000e+00_chi=15_09_May_2021_16:53',
%'Ising2D_g=2.5000e+00_chi=12_09_May_2021_17:16',
%'Ising2D_g=2.5000e+00_chi=12_09_May_2021_17:29';
%'Ising2D_g=2.5000e+00_chi=12_09_May_2021_17:36';
%'Ising2D_g=2.5000e+00_chi=8_09_May_2021_17:55';
%'Ising2D_g=2.5000e+00_chi=8_09_May_2021_18:23';
%'Ising2D_g=2.5000e+00_chi=8_09_May_2021_18:28';
%'Ising2D_g=2.5000e+00_chi=8_09_May_2021_18:33',
%'Ising2D_g=2.5000e+00_chi=12_10_May_2021_16:35',
%'Ising2D_g=2.5000e+00_chi=12_10_May_2021_17:25',
%'Ising2D_g=2.5000e+00_chi=12_10_May_2021_18:31',
%    'Ising2D_g=2.5000e+00_chi=12_10_May_2021_18:31',
%'Ising2D_T=5.0000e-01_chi=8_12_May_2021_14:55';
%'Ising2D_T=5.0000e-01_chi=8_12_May_2021_15:29';
%'Ising2D_T=5.0000e-01_chi=8_12_May_2021_16:29';
%'Ising2D_T=5.0000e-01_chi=8_12_May_2021_17:17';
%'Ising2D_T=5.0000e-01_chi=8_12_May_2021_17:53';
%'Ising2D_T=5.0000e-01_chi=16_13_May_2021_15:13'; %
%'Ising2D_T=1.0000e-01_chi=8_13_May_2021_15:52';%L=3, 20
%'Ising2D_T=1.0000e-01_chi=8_13_May_2021_17:35';%L=2
%'Ising2D_T=1.0000e-01_chi=16_13_May_2021_17:50'; % L=3, Max = 30 (tbond = 21+30)
%'Ising2D_T=1.0000e-01_chi=16_13_May_2021_19:17';%L=3, 30, chi16
%'Ising2D_T=1.0000e-01_chi=16_13_May_2021_17:59';
    }, {
%'Ising2D_g=2.5000e+00_chi=12_14_May_2021_17:07',
%'Ising2D_g=2.5000e+00_chi=12_14_May_2021_17:07',

%'Ising2D_g=2.5000e+00_chi=4_15_May_2021_16:47';
%'Ising2D_g=2.5000e+00_chi=8_15_May_2021_16:47';
%'Ising2D_g=2.5000e+00_chi=12_15_May_2021_16:47';
%'Ising2D_g=2.5000e+00_chi=16_15_May_2021_17:14';
%'Ising2D_g=2.5000e+00_chi=32_15_May_2021_18:46';

%'Ising2D_g=2.5000e+00_chi=8_16_May_2021_09:13';
%'Ising2D_g=2.5000e+00_chi=4_16_May_2021_09:10';
%'Ising2D_g=2.5000e+00_chi=4_16_May_2021_09:43';
%'Ising2D_g=2.5000e+00_chi=8_16_May_2021_09:43';

%'Ising2D_g=2.5000e+00_chi=4_16_May_2021_11:25';
%'Ising2D_g=2.5000e+00_chi=8_16_May_2021_11:25';
    'Ising2D_g=2.5000e+00_chi=8_16_May_2021_16:58';
    'Ising2D_g=2.5000e+00_chi=12_16_May_2021_16:58';
    }};

T_c_arr = [2.75, 1.2737];

%1.27375%1.2737(2)

ln_lopts = struct('Display', 0, 'maxit', 1);

filterPoints = 0;

skip = [1, 0, 1, 1];

close all

plot_opts.marker_size = 2;

for j = 1:numel(names)

    T_c = T_c_arr(j);

    if skip(j) ~= 1

        for i = 1:numel(names{j})

            data = fetch_matfiles(names{j}{i}, struct);
            data = filter_ising_results(data, struct('tol', 1e-10));

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

            marek_arr = real(data.eps_i(:, 2) - data.eps_i(:, 1));

            plot_opts.colour = marek_arr;
            chi = data.chi(1);

            %T_c = 1.2737;  %https://journals.aps.org/prb/pdf/10.1103/PhysRevB.93.155157

            plot_m_vs_t(data, chi, marek_arr, i, j, T_c, plot_opts)

            plot_m_vs_t_marek(data, chi, marek_arr, i, j, T_c, plot_opts)
            plot_xi_marek(data, chi, marek_arr, i, j, T_c, plot_opts)
            plot_S_marek(data, chi, marek_arr, i, j, T_c, plot_opts)

        end

    end
end

function plot_m_vs_t(data, chi, marek_arr, i, j, T_c, plot_opts)

    figure(j);

    y_arr = data.m;
    x_arr = data.(data.free_var);

    if i == 1
        xline(T_c, 'DisplayName', "$ T_c $")

        title(sprintf("2D transverse Ising, %s =%.4f", data.fixed_var, data.fixed_val));

        switch data.free_var
            case 'T'
                xlabel("$\frac{k T}{J}$", "Interpreter", "Latex", 'FontSize', 12);
            case 'g'
                xlabel("$\Gamma$", "Interpreter", "Latex", 'FontSize', 12);
        end

        ylabel("$\left < m \right >$", "Interpreter", "Latex");

        legend('Location', 'southwest', "Interpreter", "Latex", 'FontSize', 12)
        %colorbar
        %set(gca, 'colorscale', 'log')
        %colormap('gray')
    end

    hold on

    plot(x_arr, y_arr, '*', 'MarkerSize', plot_opts.marker_size, 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    %scatter(x_arr, y_arr, 10, plot_opts.colour, 'filled', 'DisplayName', sprintf("$ \\chi  = %d $", chi))

    hold off

end

function plot_xi_marek(data, chi, marek_arr, i, j, T_c, plot_opts)

    nu = 1;
    omega = 1;
    c = -0.5;
    d = 0.1;
    phi = 1;

    %y_arr =  (1./data.inv_corr_length .* marek_arr.^(1/nu))./( 1+c*marek_arr.^(omega)  )   ;
    %x_arr = (data.T - T_c) .* (marek_arr.^(-1/1)) + d*marek_arr.^phi/nu;

    y_arr = (1 ./ real(data.inv_corr_length) .* marek_arr.^(1 / nu));
    x_arr = (data.T - T_c) .* (marek_arr.^(-1/1));

    figure(j + 600);

    if i == 1

        title(sprintf("2D transverse Ising, %s =%.4f", data.fixed_var, data.fixed_val));

        switch data.free_var
            case 'T'
                xlabel("$\frac{k T}{J}$", "Interpreter", "Latex", 'FontSize', 12);
            case 'g'
                xlabel("$\Gamma$", "Interpreter", "Latex", 'FontSize', 12);
        end

        legend('Location', 'northwest', "Interpreter", "Latex", 'FontSize', 12)
        colorbar
        set(gca, 'colorscale', 'log')
        %colormap('gray')
    end

    hold on

    plot(x_arr, y_arr, '*', 'MarkerSize', plot_opts.marker_size, 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    %scatter(x_arr, y_arr, 10, plot_opts.colour, 'filled', 'DisplayName', sprintf("$ \\chi  = %d $", chi))

    hold off

end

function plot_S_marek(data, chi, marek_arr, i, j, T_c, plot_opts)
    c = 1/2;

    nu = 1;
    omega = 1.5;
    e = 0;
    d = 0;
    phi = 1;

    y_arr = log(exp(6 * data.S / c) .* marek_arr.^(1 / nu));
    x_arr = (data.T - T_c) .* (marek_arr.^(-1/1)) + d * marek_arr.^phi / nu;

    %y_arr = 6*data.S/c +log(marek_arr) ;
    %x_arr =(data.T- T_c).* (marek_arr.^(-1/1));

    figure(j + 900);

    if i == 1
        title(sprintf("2D transverse Ising, %s =%.4f", data.fixed_var, data.fixed_val));

        switch data.free_var
            case 'T'
                xlabel("$\frac{k T}{J}$", "Interpreter", "Latex", 'FontSize', 12);
            case 'g'
                xlabel("$\Gamma$", "Interpreter", "Latex", 'FontSize', 12);
        end

        legend('Location', 'northwest', "Interpreter", "Latex", 'FontSize', 12)
        colorbar
        set(gca, 'colorscale', 'log')
        %colormap('gray')
    end

    hold on

    plot(x_arr, y_arr, '*', 'MarkerSize', plot_opts.marker_size, 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    %scatter(x_arr, y_arr, 10, plot_opts.colour, 'filled', 'DisplayName', sprintf("$ \\chi  = %d $", chi))

end

function plot_m_vs_t_marek(data, chi, marek_arr, i, j, T_c, plot_opts)

    nu = 1;
    c = -0.3;
    omega = 1.4;
    phi = 1;
    d = -0.1;

    yarr = data.m .* (marek_arr.^(-1/8)) ./ (1);
    xarr = (data.T - T_c) / T_c .* (marek_arr.^(-1 / nu)) + d * marek_arr.^(phi / nu);

    figure(j + 400);

    if i == 1
        title(sprintf("2D transverse Ising, %s =%.4f", data.fixed_var, data.fixed_val));

        switch data.free_var
            case 'T'
                xlabel("$\frac{k T}{J}$", "Interpreter", "Latex", 'FontSize', 12);
            case 'g'
                xlabel("$\Gamma$", "Interpreter", "Latex", 'FontSize', 12);
        end

        legend('Location', 'southwest', "Interpreter", "Latex", 'FontSize', 12)

        %set(gca, 'ColorScale', 'log')

        colorbar
        set(gca, 'colorscale', 'log')
        %colormap('gray')
    end

    hold on

    plot(xarr, yarr, '*', 'MarkerSize', plot_opts.marker_size, 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    %scatter(xarr, yarr, 10, plot_opts.colour, 'filled', 'DisplayName', sprintf("$ \\chi  = %d $", chi))

    hold off

end
