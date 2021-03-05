fold = mfilename('fullpath');
pathparts = strsplit(fold, '/');

pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
fold2 = strjoin(pathparts, '/');

names = {
    'Ising2D_g=2.5000e+00_chi=2_05_March_2021_11:28.mat';
    'Ising2D_g=2.5000e+00_chi=3_05_March_2021_11:28.mat';
    'Ising2D_g=2.5000e+00_chi=4_05_March_2021_11:28.mat'};

for i = 1:numel(names)
    load(sprintf("%s/%s", fold2, names{i}));

    plot_Ising2D(m_arr, T_arr, J, g, chi)
end

function plot_Ising2D(m_arr, T_arr, J, g, chi)

    %print(f);

    figure(1);
    %loglog(  beta_arr,err_arr );
    hold on
    plot(T_arr, m_arr, '*');

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

    ylim([0, 1]);

    title(sprintf("2D transverse Ising, g=%.4f, chi=%d ", g, chi));
    xlabel("$\frac{k T}{J}$", "Interpreter", "Latex");
    ylabel("$\left < m \right >$", "Interpreter", "Latex");

    %lgd = legend('simul', 'onsager', 'fit');
    %lgd.Location = 'southwest';

    hold off

    drawnow;

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
