chi = 20;

fold = mfilename('fullpath');
pathparts = strsplit(fold, '/');

pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
fold2 = strjoin(pathparts, '/');

dt = datestr(now, 'dd_mmmm_yyyy_HH:MM');

chi_arr = [5, 10, 15, 20, 25, 30, 35, 40];
g = 2.5;

mbound = [0.6, 0.1];

calc_ising_2d(1.21, 1, g, chi_arr, 0.01, 0.01, mbound, fold2, dt, 0);

%[f,gof] = fit_data(1,1e-7,m_arr, T_arr);
%plot_onsager(m_arr, T_arr, 1,g,chi)

function calc_ising_2d(T, J, g, chi_arr, aim_dx, aim_dy, mbound, fold2, dt, onsager)

    if nargin < 6
        onsager = 0;
    end

    opts = [];
    opts.testing = 0;
    opts.visualise = 0;
    opts.double = 0;

    parfor t = 1:numel(chi_arr)
        chi = chi_arr(t);

        name = sprintf("%s/Ising2D_g=%.4e_chi=%d_%s.mat", fold2, g, chi, dt);
        T0 = T;

        d = 2;

        handle = @make_PEPO_2D_A;

        %hamiltonian setup
        S_x = [0, 1; 1, 0];
        S_y = [0, -1i; 1i, 0];
        S_z = [1, 0; 0, -1];
        I_tensor = eye(2);

        %J=1;
        %g=0.01;

        H_1_tensor = -J * g * S_x;
        H_2_tensor = -J * (reshape(ncon({S_z, S_z}, {[-1, -3], [-2, -4]}), [d, d, d, d]));

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

        m_0 = [mbound(1), mbound(1), mbound(1)];
        T_0 = [T0 - 2 * aim_dy, T0 - 1 * aim_dy, T0];

        scale_factor_dy = aim_dx / aim_dy;

        deadcounter = 1;

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
                mag = m_onsager(T0, J) / t;
                corr_len = 1;
                marek = 1;

            else

                pepo = PEPO(d, -beta * H_1_tensor, -beta * H_2_tensor, 5, handle, opts);

                [mag, inv_corr_length, marek] = PEPO_get_expectation(pepo, S_z, chi);

                corr_len = 1 / inv_corr_length;

            end

            m = abs(mag);

            fprintf("%d:%d T %.4e mag:%.4e xi:%.4e marek gap:%.4f \n", chi, i, T0, m, corr_len, marek);

            m_0(1:2) = m_0(2:3);
            m_0(3) = m;

            T_0(1:2) = T_0(2:3);
            T_0(3) = T0;

            m_arr(i) = m;
            T_arr(i) = T0;
            corr_arr(i) = corr_len;
            marek_arr(i) = marek;

            saveboy(name, 'T_arr', 'm_arr', 'corr_arr', 'marek_arr', 'chi', 'J', 'g', T_arr, m_arr, corr_arr, marek_arr, chi, J, g);

            if m < mbound(2)
                deadcounter = deadcounter -1;

                if deadcounter == 0
                    fprintf('%d finished', chi)
                    break;
                end

            end

        end

        %strip trailing places in saved array

        T_arr = T_arr(1:i);
        m_arr = m_arr(1:i);
        corr_arr = corr_arr(1:i);
        marek_arr = marek_arr(1:i);

        saveboy(name, 'T_arr', 'm_arr', 'corr_arr', 'marek_arr', 'chi', 'J', 'g', T_arr, m_arr, corr_arr, marek_arr, chi, J, g);

        fprintf("")
    end

end

function m = m_onsager(T, J)
    T_c = 2 * J / (log(1 + sqrt(2)));
    m = T;
    mask = T < T_c;
    m(mask) = (1 - sinh((2 * J) ./ T(mask)).^(-4)).^(1/8);
    m(~mask) = 0;
end
