fold = mfilename('fullpath');
pathparts = strsplit(fold, '/');

pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
fold2 = strjoin(pathparts, '/');

names = {
    'Ising2D_g=2.5000e+00_chi=15_06_March_2021_12:10.mat';
    'Ising2D_g=2.5000e+00_chi=20_06_March_2021_12:10.mat';
    'Ising2D_g=2.5000e+00_chi=25_06_March_2021_12:10.mat';
    'Ising2D_g=2.5000e+00_chi=30_06_March_2021_12:10.mat';
    'Ising2D_g=2.5000e+00_chi=35_06_March_2021_12:10.mat';
    'Ising2D_g=2.5000e+00_chi=40_06_March_2021_12:10.mat';
    };

for i = 1:numel(names)
    load(sprintf("%s/%s", fold2, names{i}));

    plot_Ising2D(m_arr, T_arr, J, 2.5, chi, i)
end

function plot_Ising2D(m_arr, T_arr, J, g, chi, i)

    %print(f);

    mask = T_arr ~= 0;
    m_arr = m_arr(mask);
    T_arr = T_arr(mask);

    figure(1);
    %loglog(  beta_arr,err_arr );
    hold on
    plot(T_arr, m_arr, '*-', 'DisplayName', sprintf("chi = %d", chi));

    %     T_min = min(T_arr);
    %     T_max = max(T_arr);
    %     T_onsager = T_min:0.001:T_max;

    %plot(T_onsager, m_onsager(T_onsager, J));

    %[f, gof] = fit_data(0.5, 1e-3, m_arr, T_arr);

    %fitfun = @(t) f.a * ((f.Tc - t) ./ f.Tc).^f.beta;
    %[f,gof] = fit_data(0.3,1e-8,   m_onsager(T_onsager,J),T_onsager  );
    %fitfun = @(t) ( (f.Tc-t)./f.Tc).^f.beta ;

    %fprintf('dbeta %.4e dTc %.4f\n', abs(f.beta - 1/8), abs(f.Tc - 2 * J / (log(1 + sqrt(2)))));

    %plot(T_arr, fitfun(T_arr));

    if i == 1
        ylim([0.25, 0.6]);
        xlim([1.23, 1.285])

        xline(1.27376)

        title(sprintf("2D transverse Ising, g=%.4f, chi=%d ", g, chi));
        xlabel("$\frac{k T}{J}$", "Interpreter", "Latex");
        ylabel("$\left < m \right >$", "Interpreter", "Latex");

        legend('Location', 'southwest')
    end
    %lgd = legend('simul', 'onsager', 'fit');
    %lgd.Location = 'southwest';

    hold off

    %drawnow;

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
