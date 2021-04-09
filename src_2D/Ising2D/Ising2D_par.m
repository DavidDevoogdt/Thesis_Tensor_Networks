function Ising2D_par(chi_arr, g, template_name)

    %for chi arr: round(2.^(3:0.5:7))

    maxit = 3;
    w = 0.5;

    switch g
        case 0
            Tk = 2.27;

            T_min = Tk - w;
            T_max = Tk + w;
        case 1.5
            Tk = 1.975;

            T_min = Tk - w;
            T_max = Tk + w;
        case 2.5
            Tk = 1.2735;

            T_min = Tk - w;
            T_max = Tk + w;
        case 2.9
            T_min = 0.4;
            T_max = 1.0;

        otherwise
            error('provide T bounds');
    end

    switch getenv('USER')
        case "david"
            nsammple = 4;
        otherwise
            nsammple = 30;
    end

    %parallel preferences
    cluster = parcluster('local');
    cluster.NumWorkers = nsammple;
    cluster.NumThreads = 1;

    %

    dt = datestr(now, 'dd_mmmm_yyyy_HH:MM');
    
    for i = 1:numel(chi_arr)

        %template setup
        if nargin == 3 %load from file
            [~, template] = fetch_matfiles(template_name, struct);
        else %make template
            fold = mfilename('fullpath');
            pathparts = strsplit(fold, '/');

            pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
            fold2 = strjoin(pathparts, '/');

            

            vumps_opts = [];
            vumps_opts.vumps_maxit = 1000;
            vumps_opts.tolfixed = 1e-10;

            template = [];

            chi = chi_arr(i);

            fprintf("\n");
            nn = sprintf("Ising2D_g=%.4e_chi=%d_%s", g, chi, dt);

            template.name = nn;

            template.name_prefix = sprintf("%s/%s", fold2, nn);

            template.T_bounds = [T_min, T_max];

            dir_name = sprintf("%s/", template.name_prefix);
            if ~exist(dir_name, 'dir')
                mkdir(dir_name);
            end

            %make template for simulations
            template.vumps_opts = vumps_opts;
            template.handle = @make_PEPO_2D_B;

            template.dir_name = dir_name;
            template.time = dt;

            S_z = [1, 0; 0, -1];
            template.X = S_z; %observable

            template.model_params = models('t_ising', struct('g', g));
            template.pepo_opts = struct();

            saveboy(sprintf("%s/template.mat", template.name_prefix), 'template', template);
        end

        disp(template.name_prefix)

        template.vumps_opts.chi_max = chi;
        calc_ising_2d(T_min, T_max, 0.01, 0.01, nsammple, template, maxit);
    end
end

function calc_ising_2d(Tmin, Tmax, aim_dx, aim_dy, nsammple, template, maxit)

    for i = 1:maxit

        %determine next temperatures
        if i == 1
            T0 = (Tmax - Tmin) / (nsammple - 1) * (0:nsammple - 1) + Tmin;
        else

            data = fetch_matfiles(template.name, struct);
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
