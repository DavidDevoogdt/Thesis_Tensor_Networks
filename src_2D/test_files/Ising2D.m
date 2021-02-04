[m_arr, T_arr, corr_arr, marek_arr, name] = calc_ising_2d(0.4, 1.0, 1 , 3, 0.05, 0.02, 0);

%[f,gof] = fit_data(1,1e-7,m_arr, T_arr);
plot_onsager(m_arr, T_arr, 1, 0.001)

function [m_arr, T_arr, corr_arr, marek_arr, name] = calc_ising_2d(T0, J, g, chi, aim_dx, aim_dy, onsager)

    if nargin < 5
        onsager = 0;
    end

    if onsager == 1
        name = sprintf("onsager.mat");
    else
        name = sprintf("Ising2D_g=%.4e_chi=%d.mat", g, chi);
    end

    d = 2;

    handle = @make_PEPO_2D_B;
    
    %hamiltonian setup
    S_x = [0, 1; 1, 0];
    S_y = [0, -1i; 1i, 0];
    S_z = [1, 0; 0, -1];
    I_tensor = eye(2);

    %J=1;
    %g=0.01;

    H_1_tensor = -J * g * S_x;
    H_2_tensor = -J * (reshape(ncon({S_z, S_z}, {[-1, -3], [-2, -4]}), [d, d, d, d]));

    opts.testing = 0;
    opts.visualise = 0;

    pos_map = [0, 0, 1, 1;
                1, 1, 1, 1;
                0, 0, 1, 1];

    map = create_map(pos_map);

    T_c = 2 * J / (log(1 + sqrt(2)));

    maxit = 100;

    T_arr = zeros(1, maxit);
    m_arr = zeros(1, maxit);
    corr_arr = zeros(1, maxit);
    marek_arr = zeros(1, maxit);

    %T0=1;

    %aim_dx = 0.1;
    %aim_dy = 0.02;%tends to overshoot, make small enough
    ds = max(aim_dx, aim_dy);

    m_0 = [1, 1, 1];
    T_0 = [T0 - 2 * aim_dy, T0 - 1 * aim_dy, T0];

    scale_factor_dy = aim_dx / aim_dy;

    deadcounter = 2;

    for i = 1:maxit
        p = polyfit(T_0, m_0, 2);
        dp = polyder(p);
        dmdT = @(x) polyval(dp, x);

        f = @(t) sqrt(1 + (scale_factor_dy * (dmdT(t))).^2);
        arc_len = @(dt) integral(f, T0, T0 + dt) - ds;
        dT = fzero(arc_len, [0, ds * 1.01]);

        T0 = T0 + dT;

        beta = 1 / T0;

        if onsager == 1
            mag = m_onsager(T0, J);
            corr_len = 1;
            marek = 1;

        else
            pepo = PEPO(d, -beta * H_1_tensor, -beta * H_2_tensor, 5, handle, opts);


            
            [mag, inv_corr_length, marek] = PEPO_get_expectation(pepo,S_z, chi);

            corr_len = 1 / inv_corr_length;

        end

        m = abs(mag);

        fprintf(" T %.4e mag %.4e ,Theory: %.4e  corrlen %.4e marek gap %.4f \n", T0, m, m_onsager(T0, J), corr_len, marek);

        m_0(1:2) = m_0(2:3);
        m_0(3) = m;

        T_0(1:2) = T_0(2:3);
        T_0(3) = T0;

        m_arr(i) = m;
        T_arr(i) = T0;
        corr_arr(i) = corr_len;
        marek_arr(i) = marek;

        if m < 1e-5
            deadcounter = deadcounter -1;

            if deadcounter == 0
                break;
            end

        end

    end

    T_arr = T_arr(1:i);
    m_arr = m_arr(1:i);
    corr_arr = corr_arr(1:i);
    marek_arr = marek_arr(1:i);

    save(name, 'T_arr', 'm_arr', 'corr_arr', 'marek_arr');

end

function plot_onsager(m_arr, T_arr, J, g)

    %print(f);

    figure();
    %loglog(  beta_arr,err_arr );
    plot(T_arr, m_arr, '*');

    hold on
    T_min = min(T_arr);
    T_max = max(T_arr);
    T_onsager = T_min:0.001:T_max;

    plot(T_onsager, m_onsager(T_onsager, J));

    [f, gof] = fit_data(0.7, 1e-1, m_arr, T_arr);

    fitfun = @(t) f.a * ((f.Tc - t) ./ f.Tc).^f.beta;
    %[f,gof] = fit_data(0.3,1e-8,   m_onsager(T_onsager,J),T_onsager  );
    %fitfun = @(t) ( (f.Tc-t)./f.Tc).^f.beta ;

    fprintf('dbeta %.4e dTc %.4f\n', abs(f.beta - 1/8), abs(f.Tc - 2 * J / (log(1 + sqrt(2)))));

    plot(T_arr, fitfun(T_arr));

    ylim([0, 1]);

    title(sprintf("2D transverse Ising, g=%.4f ", g));
    xlabel("$\frac{k T}{J}$", "Interpreter", "Latex");
    ylabel("$\left < m \right >$", "Interpreter", "Latex");

    legend('simul', 'onsager', 'fit')

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

function m = m_onsager(T, J)
    T_c = 2 * J / (log(1 + sqrt(2)));
    m = T;
    mask = T < T_c;
    m(mask) = (1 - sinh((2 * J) ./ T(mask)).^(-4)).^(1/8);
    m(~mask) = 0;
end
