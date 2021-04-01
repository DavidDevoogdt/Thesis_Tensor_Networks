%notes: http://alps.comp-phys.org/mediawiki/index.php/Main_Page
% https://arxiv.org/pdf/1406.2973.pdf

%sizes: 20*1.2.^(1:9)

% names = {{
% %'Ising2D_g=2.5000e+00_chi=24_15_March_2021_15:41';
%     'Ising2D_g=2.5000e+00_chi=15_20_March_2021_17:45';
%     'Ising2D_g=2.5000e+00_chi=29_15_March_2021_15:41';
%     'Ising2D_g=2.5000e+00_chi=20_20_March_2021_12:48';
%     'Ising2D_g=2.5000e+00_chi=35_15_March_2021_15:41';
%     'Ising2D_g=2.5000e+00_chi=41_15_March_2021_15:41';
%     'Ising2D_g=2.5000e+00_chi=50_15_March_2021_16:20';
%     'Ising2D_g=2.5000e+00_chi=60_15_March_2021_16:54';
%     'Ising2D_g=2.5000e+00_chi=73_15_March_2021_19:35';
%     'Ising2D_g=2.5000e+00_chi=86_15_March_2021_19:46';
%     }};

names = {{
%     'Ising2D_g=2.5000e+00_chi=20_22_March_2021_10:11_rep';
%     'Ising2D_g=2.5000e+00_chi=25_22_March_2021_10:11_rep';
%     'Ising2D_g=2.5000e+00_chi=30_22_March_2021_10:11_rep';
%     'Ising2D_g=2.5000e+00_chi=40_22_March_2021_10:46_rep';
%     'Ising2D_g=2.5000e+00_chi=60_22_March_2021_11:49_rep';
%     'Ising2D_g=2.5000e+00_chi=65_22_March_2021_13:52_rep';
%     'Ising2D_g=2.5000e+00_chi=70_22_March_2021_13:51_rep';
%     }};
%     },{
%'Ising2D_g=2.5000e+00_chi=20_22_March_2021_10:11';
%'Ising2D_g=2.5000e+00_chi=25_22_March_2021_10:11';
%'Ising2D_g=2.5000e+00_chi=30_22_March_2021_10:11';
    'Ising2D_g=2.5000e+00_chi=40_22_March_2021_10:46';
    'Ising2D_g=2.5000e+00_chi=60_22_March_2021_11:49';
    'Ising2D_g=2.5000e+00_chi=65_22_March_2021_13:52';
    'Ising2D_g=2.5000e+00_chi=70_22_March_2021_13:51';
    }};

ln_lopts = struct('Display', 0, 'maxit', 1);

filterPoints = 1;

close all

for j = 1:numel(names)
    for i = 1:numel(names{j})

        data = fetch_matfiles(names{j}{i}, struct);
        data = filter_ising_results(data, struct);

        marek_arr = data.marek;
        corr_arr = 1 ./ data.inv_corr_length;
        m_arr = data.m;
        T_arr = data.T;
        J = data.J;
        g = data.g;
        chi = data.chi;

        T_c = 1.2736;
        %T_c = 1.2737;  %https://journals.aps.org/prb/pdf/10.1103/PhysRevB.93.155157

        plot_m_vs_t(m_arr, T_arr, J, g, chi, i, j, T_c)
        plot_m_vs_t_marek(marek_arr, corr_arr, m_arr, T_arr, J, g, chi, i, j, T_c)

        plot_xi(marek_arr, corr_arr, m_arr, T_arr, J, g, chi, i, j, T_c)
        plot_xi_marek(marek_arr, corr_arr, m_arr, T_arr, J, g, chi, i, j, T_c)
        %crit_plot(marek_arr, corr_arr, m_arr, T_arr, J, g, chi, i, j, T_c)
        %collapse(corr_arr, m_arr, T_arr, J, g, chi, i, j)

        %delta_scatter(marek_arr, m_arr, T_arr, J, g, chi, i, j, T_c);

    end
end

% function delta_scatter(marek_arr, m_arr, T_arr, J, g, chi, i, j, T_c)

%     figure(j + 700);

%     % if i == 1
%     %     %     xline(1.27376, 'DisplayName', "$ T_c $")

%     %     %     title(sprintf("2D transverse Ising, g=%.4f", g));
%     %     %     xlabel("$\frac{k T}{J}$", "Interpreter", "Latex", 'FontSize', 12);
%     %     %     ylabel("$\left < m \right >$", "Interpreter", "Latex");

%     %     %     legend('Location', 'southwest', "Interpreter", "Latex", 'FontSize', 12)
%     %     colorbar
%     %     set(gca, 'ColorScale', 'log')
%     % end

%     hold on
%     plot3(T_arr, m_arr, marek_arr);

%     hold off

% end

function plot_m_vs_t(m_arr, T_arr, J, g, chi, i, j, T_c)

    figure(j);

    if i == 1
        xline(1.27376, 'DisplayName', "$ T_c $")

        title(sprintf("2D transverse Ising, g=%.4f", g));
        xlabel("$\frac{k T}{J}$", "Interpreter", "Latex", 'FontSize', 12);
        ylabel("$\left < m \right >$", "Interpreter", "Latex");

        legend('Location', 'southwest', "Interpreter", "Latex", 'FontSize', 12)
    end

    hold on
    plot(T_arr, m_arr, '*-', 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    hold off

end

function plot_xi(marek_arr, corr_arr, m_arr, T_arr, J, g, chi, i, j, T_c)

    y_arr = corr_arr;
    x_arr = T_arr;

    figure(j + 100);

    if i == 1

        title(sprintf("2D transverse Ising, g=%.4f", g));
        xlabel("$T$", "Interpreter", "Latex", 'FontSize', 12);
        ylabel("$ \xi  $", "Interpreter", "Latex");

        legend('Location', 'northwest', "Interpreter", "Latex", 'FontSize', 12)
    end

    hold on
    plot(x_arr, real(y_arr), '*-', 'DisplayName', sprintf("$ \\chi  = %d $", chi));
    hold off

end

function plot_xi_marek(marek_arr, corr_arr, m_arr, T_arr, J, g, chi, i, j, T_c)

    y_arr = corr_arr .* marek_arr;
    x_arr = (T_arr - T_c) .* (marek_arr.^(-1/1));

    figure(j + 600);

    if i == 1

        title(sprintf("2D transverse Ising, g=%.4f", g));
        xlabel("$(T-T_c)\delta^{-1/\nu}$", "Interpreter", "Latex", 'FontSize', 12);
        ylabel("$ \xi \delta  $", "Interpreter", "Latex");

        %legend('Location', 'northwest', "Interpreter", "Latex", 'FontSize', 12)

        colorbar
        %set(gca, 'colorscale', 'log')
        %colormap('gray')
        % cb = colorbar();
        % cb.Ruler.Scale = 'log';
        % cb.Ruler.MinorTick = 'on';
    end

    hold on
    %plot(x_arr, real(y_arr), '*', 'DisplayName', sprintf("$ \\chi  = %d $", chi));
    scatter(x_arr, real(y_arr), 10, marek_arr, 'filled')

    hold off

end

function crit_plot(marek_arr, corr_arr, m_arr, T_arr, J, g, chi, i, j, T_c)

    %print(f);

    yarr = m_arr .* marek_arr.^(-1/8);
    xarr = T_arr;

    %yarr = m_arr .* marek_arr.^(-1/8);
    %xarr = (T_arr-T_c).*marek_arr.^(-1);

    figure(j + 200);
    %loglog(  beta_arr,err_arr );

    if i == 1
        ylim([0, 1]);
        xlim([T_c - 0.01, T_c + 0.01])

        xline(1.27376, 'DisplayName', "$ T_c $")

        title(sprintf("2D transverse Ising, g=%.4f", g));
        xlabel("$\frac{k T}{J}$", "Interpreter", "Latex", 'FontSize', 12);
        ylabel("$ m \xi ^{\frac{\beta}{\nu}}  $", "Interpreter", "Latex");

        legend('Location', 'northwest', "Interpreter", "Latex", 'FontSize', 12)
    end

    hold on
    plot(xarr, yarr, '*', 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    idxmax = find(yarr < 0.45);
    Tmax = T_arr(idxmax(1));

    Tmin = T_arr(1);

    fict_T = Tmin:1e-4:Tmax;
    vq = interp1(T_arr(1:idxmax), yarr(1:idxmax), fict_T, 'spline');
    plot(fict_T, vq, '-', 'DisplayName', sprintf("$ interp \\chi  = %d $", chi));

    %lgd = legend('simul', 'onsager', 'fit');
    %lgd.Location = 'southwest';

    hold off

    %drawnow;

end

function collapse(corr_arr, m_arr, T_arr, J, g, chi, i, j, T_c)

    %print(f);

    yarr = m_arr .* corr_arr.^(1/8);
    xarr = (T_arr - T_c) ./ T_c .* corr_arr;

    figure(j + 300);
    %loglog(  beta_arr,err_arr );
    hold on
    plot(xarr, yarr, '*', 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    if i == 1
        ylim([0, 1]);
        xlim([-0.5, 0.6])

        xline(1.27376, 'DisplayName', "$ T_c $")

        title(sprintf("2D transverse Ising, g=%.4f", g));
        xlabel("$t \xi^{\frac{1}{\nu}}$", "Interpreter", "Latex", 'FontSize', 12);
        ylabel("$ m \xi ^{\frac{\beta}{\nu}}  $", "Interpreter", "Latex");

        legend('Location', 'southwest', "Interpreter", "Latex", 'FontSize', 12)
    end
    %lgd = legend('simul', 'onsager', 'fit');
    %lgd.Location = 'southwest';

    hold off

    %drawnow;

end

function plot_m_vs_t_marek(marek_arr, corr_arr, m_arr, T_arr, J, g, chi, i, j, T_c)

    yarr = m_arr .* (marek_arr.^(-1/8));
    xarr = (T_arr - T_c) / T_c .* (marek_arr.^(-1));

    figure(j + 400);

    if i == 1
        title(sprintf("2D transverse Ising, g=%.4f ", g));
        xlabel("$ (T-T_C) \delta^{-1 / \nu} $", "Interpreter", "Latex", 'FontSize', 12);
        ylabel("$  \left< m \right>  \delta ^{ -\beta / \nu}  $", "Interpreter", "Latex");

        %legend('Location', 'southwest', "Interpreter", "Latex", 'FontSize', 12)

        %set(gca, 'ColorScale', 'log')

        % c1 = min(min(marek_arr));
        % c2 = max(max(marek_arr));
        % % set limits for the caxis
        % caxis([log10(c1) log10(c2)]);
        % % preallocate Ticks and TickLabels
        % num_of_ticks = 5;
        % Ticks = zeros(1, num_of_ticks);
        % TickLabels = zeros(1, num_of_ticks);
        % % distribute Ticks and TickLabels
        % for n = 1:1:num_of_ticks

        %     Ticks(n) = log10(round(c2) / num_of_ticks * n);
        %     TickLabels(n) = round(c2) / num_of_ticks * n;
        % end
        % % set Ticks and TickLabels
        % colorbar('Ticks', Ticks, 'TickLabels', TickLabels)
        colorbar
        set(gca, 'colorscale', 'log')
        %colormap('gray')
    end

    hold on

    %plot3(xarr, yarr, marek_arr, '*', 'DisplayName', sprintf("$ \\chi  = %d $", chi));

    scatter(xarr, yarr, 10, marek_arr, 'filled')

    hold off

end

%%

function [f, gof_info] = fit_data(m_max, m_min, m_arr, T_arr)

    m_mask = (m_arr > m_min) & (m_arr < m_max);

    T_data = reshape(T_arr(m_mask), [], 1);
    m_data = reshape(m_arr(m_mask), [], 1);

    [f, gof_info] = fit(m_data, T_data, ...
        'Tc*(1- (x/a)^(1/beta) )', ...
        'StartPoint', [2.26, 1, 1/8]);

    %      [f, gof_info] = fit(m_data,T_data,...
    %         'Tc*(1- (x)^(1/beta) )',...
        %         'StartPoint', [2.26,  1/8]);

end
