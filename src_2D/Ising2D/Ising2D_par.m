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
            calc_ising_2d(1.26, 1.29, 1, g, chi_arr, 0.01, 0.01, fold2, dt, 0, nsammple, mbound, vumps_maxit);
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

        template.name_prefix = sprintf("%s/%s", fold2, nn);

        dir_name = sprintf("%s/", template.name_prefix);
        if ~exist(dir_name, 'dir')
            mkdir(dir_name);
        end

        maxit = 30;

        %make template for simulations

        vumps_opts.chi_max = chi;
        vumps_opts.vumps_maxit = vumps_maxit;
        vumps_opts.tolfixed = 1e-10;

        template.vumps_opts = vumps_opts;
        template.handle = @make_PEPO_2D_A;

        template.dir_name = dir_name;
        template.time = dt;

        S_z = [1, 0; 0, -1];
        template.X = S_z; %observable

        template.model_params = models('t_ising', struct('g', 2.5));
        template.pepo_opts = struct();

        saveboy(sprintf("%s/template_%s.mat", template.name_prefix, template.name_suffix), 'template', template);

        for i = 1:maxit

            %determine next temperatures
            if i == 1
                T0 = (Tmax - Tmin) / (nsammple - 1) * (0:nsammple - 1) + Tmin;
            else

                data = fetch_matfiles(nn, struct);
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

            parfor j = 1:numel(T0)
                save_vars = [];
                save_vars.fname = sprintf('%2d:%2d', i, j);

                Ising2D_core(save_vars, template, T0(j));
            end
        end
    end

end
