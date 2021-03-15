function Ising2D_par(chi, name)

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');

    pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
    fold2 = strjoin(pathparts, '/');

    dt = datestr(now, 'dd_mmmm_yyyy_HH:MM');

    %chi_arr = [5, 10, 15];

    mbound = [0.1, 1]; %no longer used

    g = 2.5;
    vumps_maxit = 1500;

    switch getenv('USER')
        case "david"
            chi_arr = [chi];
            nsammple = 15;
            calc_ising_2d(1.26, 1.28, 1, g, chi_arr, 0.01, 0.01, fold2, dt, 0, nsammple, mbound, vumps_maxit);
            %calc_ising_2d(2.1, 2.3, 1, g, chi_arr, 0.01, 0.01, fold2, dt, 1, nsammple, mbound, vumps_maxit);
        otherwise
            nsammple = 15;
            chi_arr = [chi];
            calc_ising_2d(1.26, 1.28, 1, g, chi_arr, 0.01, 0.02, fold2, dt, 0, nsammple, mbound, vumps_maxit);
    end
end

%

%[f,gof] = fit_data(1,1e-7,m_arr, T_arr);
%plot_onsager(m_arr, T_arr, 1,g,chi)

function calc_ising_2d(Tmin, Tmax, J, g, chi_arr, aim_dx, aim_dy, fold2, dt, onsager, nsammple, mbound, vumps_maxit)

    if nargin < 6
        onsager = 0;
    end

    datet = dt;

    cluster = parcluster('local');
    %cluster = parpool('threads')

    cluster.NumWorkers = nsammple;
    cluster.NumThreads = 1;

    for t = 1:numel(chi_arr)
        chi = chi_arr(t);

        cc = sprintf("%s/Ising2D_g=%.4e_chi=%d_%s", fold2, g, chi, datet);
        fprintf("%s.mat", cc);

        maxit = 30;

        T_arr = zeros(maxit, nsammple);
        m_arr = zeros(maxit, nsammple);
        corr_arr = zeros(maxit, nsammple);
        marek_arr = zeros(maxit, nsammple);
        ctr_arr = zeros(maxit, nsammple);
        vumps_err_arr = zeros(maxit, nsammple);
        pepo_arr = cell(maxit, nsammple);

        %T0=1;

        for i = 1:maxit

            if i == 1
                T0 = (Tmax - Tmin) / (nsammple - 1) * (0:nsammple - 1) + Tmin;
            else

                [T_arr_3, m_arr_3, corr_arr_3, marek_arr_3, vumps_err_arr_3, ~, ~, ~] = fetch_matfiles(cc);
                [T_arr_2, m_arr_2, ~, ~] = filter_ising_results(T_arr_3, m_arr_3, marek_arr_3, corr_arr_3, vumps_err_arr_3, 1);

                %idx2 = find(m_arr_2 < mbound(1));

                dT_arr = diff(T_arr_2);
                dm_arr = diff(m_arr_2);

                %dT_arr = dT_arr(1:idx2 - 1);
                %dm_arr = dm_arr(1:idx2 - 1);

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

            T_arr(i, :) = T0;

            [m_arr, T_arr, corr_arr, marek_arr, ctr_arr, vumps_err_arr,pepo_arr] = Ising2D_core(cc, J, g, onsager, chi, i, m_arr, T_arr, corr_arr, marek_arr, ctr_arr, vumps_err_arr, vumps_maxit,pepo_arr);

        end

        fprintf("")
    end

end
