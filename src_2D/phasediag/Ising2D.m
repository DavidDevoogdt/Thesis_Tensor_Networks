function Ising2D_par(chi_arr, fixed_val, fixed_var, opts, pepo_opts)
    %routine to calculate a phase diagram for transversal ising model in 2D
    %samples in batch points along m-T / m-g graph such that arc length between points is small
    %example: Ising2D_par(8, 2.5, 'g', struct('testing',1,'unit_cell',1,'par',0 ))
    %Ising2D_par(8, 0.7, 'T', struct('testing',0,'unit_cell',1,'par',1,) , struct('order',6,'max_bond_dim',20 ,'do_loops',1 ));
    %Ising2D_par(6, [], [],
    %struct('template_name','TIM_g=2.5_order_5_chi=6_sym=1_02_August_2021_15:36'))
    %(continue with exsisting set of data points)
    %for PEPO_opts, see PEPO.m

    if nargin < 5
        pepo_opts = struct(...
            'testing', 0, ...
            'order', 5);
    end

    %parse opts
    p = inputParser;
    addParameter(p, 'template_name', [])
    addParameter(p, 'par', 1)
    addParameter(p, 'sym', 1)
    addParameter(p, 'unit_cell', 1)
    addParameter(p, 'testing', 0)
    addParameter(p, 'x_bounds', [])
    addParameter(p, 'nsample', -1)
    addParameter(p, 'npoints', 70)
    parse(p, opts)

    %parallel preferences
    if p.Results.nsample == -1
        switch getenv('USER')
            case "david"
                nsample = 4;
            otherwise
                nsample = 20;
        end
    else
        nsample = p.results.nsmaple;
    end

    if p.Results.par == 0
        nsample = 1;
    end

    %cluster = parcluster('local');
    %cluster.NumWorkers = nsample;
    %cluster.NumThreads = 1;

    if isempty(gcp('nocreate'))
        parpool(nsample);
    end

    for i = 1:numel(chi_arr)

        %template setup
        if ~isempty(p.Results.template_name) %load from file
            [~, template] = sampling_fetch(p.Results.template_name, struct);
            chi = chi_arr(i);
            first = 0;
        else %make template
            first = 1;
            template = [];
            %simulation stuff
            chi = chi_arr(i);
            if isempty(p.Results.x_bounds)
                w = 0.1;
                switch fixed_var
                    case 'g'
                        switch fixed_val
                            case 0
                                Tk = 2.27;

                                x_min = Tk - w;
                                x_max = Tk + w;
                            case 1.5
                                Tk = 1.975;

                                x_min = Tk - w;
                                x_max = Tk + w;
                            case 2.5
                                Tk = 1.2735;

                                x_min = Tk - w;
                                x_max = Tk + w;
                            case 2.9
                                x_min = 0.5;
                                x_max = 0.75;

                            otherwise
                                error('provide T bounds');
                        end
                    case 'T'
                        switch fixed_val
                            case 1.2737
                                x_min = 2.3;
                                x_max = 2.7;
                            case 1
                                x_min = 2.3;
                                x_max = 2.7;
                            case 0.7
                                x_min = 2.7;
                                x_max = 2.9;
                            case 0.5
                                x_min = 2;
                                x_max = 3.0;

                            case 0.1
                                x_min = 2.5;
                                x_max = 3.0;

                            otherwise
                                error('todo')
                        end

                    otherwise
                        error('unknown')
                end
                template.x_bounds = [x_min, x_max];
            else
                template.x_bounds = p.Results.x_bounds;
            end

            template.fixed_var = fixed_var;
            template.fixed_val = fixed_val;
            template.model_name = 't_ising';
            switch fixed_var
                case 'g'
                    template.free_var = 'T';
                case 'T'
                    template.free_var = 'g';
                otherwise
                    error('todo')
            end

            if p.Results.sym == 1
                template.handle = @make_PEPO_2D_sym;
            else
                template.handle = @make_PEPO_2D_asym;
            end

            template.pepo_opts = pepo_opts;

            %folder and nam,ing
            dt = datestr(now, 'dd_mmmm_yyyy_HH:MM');

            fold = mfilename('fullpath');
            pathparts = strsplit(fold, '/'); pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
            fold2 = strjoin(pathparts, '/');

            if ~isfield(pepo_opts, 'max_bond_dim')
                nn = sprintf("TIM_%s=%.1f_order_%d_chi=%d_sym=%d_%s", fixed_var, fixed_val, pepo_opts.order, chi, p.Results.sym, dt);

            else
                nn = sprintf("TIM_%s=%.1f_order_%d_chi=%d_trunc_%d_sym=%d_%s", fixed_var, fixed_val, pepo_opts.order, chi, pepo_opts.order.max_bond_dim, p.Results.sym, dt);

            end

            template.name = nn;
            template.name_prefix = sprintf("%s/%s", fold2, nn);

            dir_name = sprintf("%s/", template.name_prefix);
            if ~exist(dir_name, 'dir')
                mkdir(dir_name);
            end

            template.dir_name = dir_name;
            template.time = dt;

            % VUMPS
            vumps_opts = [];
            vumps_opts.vumps_maxit = 1000;
            vumps_opts.tolfixed = 1e-10;

            vumps_opts.cell_size = p.Results.unit_cell;

            if ~p.Results.testing
                vumps_opts.disp = 'None';
            else
                vumps_opts.disp = 'iter';
            end

            S_z = [1, 0; 0, -1];
            template.X = S_z; %observable
            template.vumps_opts = vumps_opts;

            %
            saveboy(sprintf("%s/template.mat", template.name_prefix), 'template', template);

        end

        disp(template.name_prefix)

        template.vumps_opts.chi_max = chi;
        sampling_workers(0.01, 0.01, nsample, template, p.Results.npoints, first, p.Results.par);
    end
end
