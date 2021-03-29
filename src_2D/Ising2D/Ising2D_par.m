function Ising2D_par(chi, name)

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');

    pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
    fold2 = strjoin(pathparts, '/');

    dt = datestr(now, 'dd_mmmm_yyyy_HH:MM');

    %chi_arr = [5, 10, 15];

    mbound = [0.1, 1]; %no longer used

    g = 2.5;
    vumps_maxit = 1300;

    switch getenv('USER')
        case "david"
            chi_arr = [chi];
            nsammple = 4;
            calc_ising_2d(1.26, 1.28, 1, g, chi_arr, 0.01, 0.01, fold2, dt, 0, nsammple, mbound, vumps_maxit);
            %calc_ising_2d(2.1, 2.3, 1, g, chi_arr, 0.01, 0.01, fold2, dt, 1, nsammple, mbound, vumps_maxit);
        otherwise
            nsammple = 15;
            chi_arr = [chi];
            calc_ising_2d(1.27, 1.275, 1, g, chi_arr, 0.01, 0.02, fold2, dt, 0, nsammple, mbound, vumps_maxit);
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

    fprintf("\n");

    for t = 1:numel(chi_arr)
        chi = chi_arr(t);

        nn = sprintf("Ising2D_g=%.4e_chi=%d_%s", g, chi, datet);
        cc = sprintf("%s/%s", fold2, nn);
        fprintf("%s.mat", cc);

        maxit = 30;

        for i = 1:maxit

            if i == 1
                T0 = (Tmax - Tmin) / (nsammple - 1) * (0:nsammple - 1) + Tmin;
            else

                data = fetch_matfiles(nn,struct);
                data = filter_ising_results(data, struct);
                
                T_arr_2 = data.T;
                m_arr_2 = data.m;

                dT_arr = diff(T_arr_2);
                dm_arr = diff(m_arr_2);

                if sum(abs(dT_arr) > aim_dx) == 0 && sum(abs(dm_arr) > aim_dy) == 0
                    break;
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
          
            template.J=J;
            template.g=g;
            template.chi=chi;
            template.vumps_maxit = vumps_maxit;
            
            
            Ising2D_core(cc, T0,template,i);
            

        end

        fprintf("")
    end

end
