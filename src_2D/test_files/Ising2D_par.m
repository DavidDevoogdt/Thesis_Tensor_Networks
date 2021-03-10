fold = mfilename('fullpath');
pathparts = strsplit(fold, '/');

pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
fold2 = strjoin(pathparts, '/');

dt = datestr(now, 'dd_mmmm_yyyy_HH:MM');

%chi_arr = [5, 10, 15];

mbound = [0.3, 1];

g = 2.5;

nsammple = 10;

switch getenv('USER')
    case "david"
        chi_arr = [15];
        calc_ising_2d(1.1, 2.8, 1, g, chi_arr, 0.1, 0.1, fold2, dt, 1, nsammple, mbound);
    otherwise
        nsammple = 20;
        chi_arr = [45, 50];
        calc_ising_2d(1.26, 1.28, 1, g, chi_arr, 0.05, 0.01, fold2, dt, 0, nsammple, mbound);
end

%

%[f,gof] = fit_data(1,1e-7,m_arr, T_arr);
%plot_onsager(m_arr, T_arr, 1,g,chi)

function calc_ising_2d(Tmin, Tmax, J, g, chi_arr, aim_dx, aim_dy, fold2, dt, onsager, nsammple, mbound)

    if nargin < 6
        onsager = 0;
    end

    datet = dt;

    T_max = inf;

    cluster = parcluster('local');
    %cluster = parpool('threads')

    cluster.NumWorkers = nsammple;
    cluster.NumThreads = 1;

    for t = 1:numel(chi_arr)
        chi = chi_arr(t);

        name = sprintf("%s/Ising2D_g=%.4e_chi=%d_%s.mat", fold2, g, chi, datet);
        disp(name);

        maxit = 100;

        T_arr = zeros(1, maxit * nsammple);
        m_arr = zeros(1, maxit * nsammple);
        corr_arr = zeros(1, maxit * nsammple);
        marek_arr = zeros(1, maxit * nsammple);

        %T0=1;

        for i = 1:maxit
            if i == 1
                T0 = (Tmax - Tmin) / (nsammple - 1) * (0:nsammple - 1) + Tmin;
            else

                T0 = [];

                mask = T_arr > 0;

                [T_arr_2, idx] = sort(T_arr(mask));
                m_arr_2 = m_arr(mask);
                m_arr_2 = m_arr_2(idx);

                idx2 = find(m_arr_2 < mbound(1));

                dT_arr = diff(T_arr_2);
                dm_arr = diff(m_arr_2);

                dT_arr = dT_arr(1:idx2 - 1);
                dm_arr = dm_arr(1:idx2 - 1);

                if sum(abs(dT_arr) > aim_dx) == 0 && sum(abs(dm_arr) > aim_dy) == 0
                    break; %goal achieved
                end

                ds = (dT_arr.^2 + ((aim_dx / aim_dy) .* dm_arr).^2).^(0.5);

                new_T = zeros(size(ds));

                for ii = 1:nsammple
                    [~, idx] = max(ds ./ (new_T + 1));
                    new_T(idx) = new_T(idx) + 1;
                end

                T0 = [];

                K = find(new_T ~= 0);
                for ii = 1:numel(K)
                    ind = K(ii);
                    nnt = new_T(ind);
                    nt = T_arr_2(ind);
                    dt = dT_arr(ind) / (nnt + 1);

                    T0 = [T0, nt + dt * (1:nnt)];
                end

            end

            [m0, corr_len0, marek0, T_max] = Ising2D_core(T0, J, g, onsager, chi, i, T_max);

            m_arr((i - 1) * nsammple + 1:i * nsammple) = m0;
            T_arr((i - 1) * nsammple + 1:i * nsammple) = T0;
            corr_arr((i - 1) * nsammple + 1:i * nsammple) = corr_len0;
            marek_arr((i - 1) * nsammple + 1:i * nsammple) = marek0;

            saveboy(name, 'T_arr', 'm_arr', 'corr_arr', 'marek_arr', 'chi', 'J', 'g', T_arr, m_arr, corr_arr, marek_arr, chi, J, g);

            %print current status
            %mask = T_arr > 0;
            %[T_arr_2, idx] = sort(T_arr(mask));
            %m_arr_2 = m_arr(mask);
            %m_arr_2 = m_arr_2(idx);
            %disp(T_arr_2);
            %disp(m_arr_2);
            %
        end

        %strip trailing places in saved array

        T_arr = T_arr(1:i * nsammple);
        m_arr = m_arr(1:i * nsammple);
        corr_arr = corr_arr(1:i * nsammple);
        marek_arr = marek_arr(1:i * nsammple);

        saveboy(name, 'T_arr', 'm_arr', 'corr_arr', 'marek_arr', 'chi', 'J', 'g', T_arr, m_arr, corr_arr, marek_arr, chi, J, g);

        fprintf("")
    end

end
