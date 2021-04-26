function Ising2D_par(chi_arr, fixed_val, fixed_var ,template_name)

    %for chi arr: round(2.^(3:0.5:7))

    maxit = 5;
    w = 0.5;

    
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
                    x_min = 0.4;
                    x_max = 1.0;

                otherwise
                    error('provide T bounds');
            end    
        case 'T'
            switch fixed_val
                case 1.2737
                    x_min = 2.3;
                    x_max = 2.7;
                
                case 0.7
                    x_min = 2.6;
                    x_max = 2.9;
                case 0.5
                    x_min = 0;
                    x_max = 3.0;
                    
                 case 0.1
                    x_min = 0;
                    x_max = 3.0;
                    
                    
                otherwise
                    error('todo')
            end
                    
        otherwise
            error('unknown')
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
        if nargin == 4 %load from file
            [~, template] = fetch_matfiles(template_name{i}, struct);
            chi = chi_arr(i);
            
            first=0;
            
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
            
            nn = sprintf("Ising2D_%s=%.4e_chi=%d_%s", fixed_var ,fixed_val , chi, dt);

            template.name = nn;
            template.name_prefix = sprintf("%s/%s", fold2, nn);

            template.x_bounds = [x_min, x_max];

            template.fixed_var =  fixed_var;
            switch fixed_var
                case 'g'
                    template.free_var = 'T';
                case 'T'
                    template.free_var = 'g';
                otherwise
                    error('todo')
            end
            
   
            template.fixed_val =  fixed_val;
            
            
            template.model_name = 't_ising';
            
            
            %template.model_params = models('t_ising', struct('g', g));
            
            dir_name = sprintf("%s/", template.name_prefix);
            if ~exist(dir_name, 'dir')
                mkdir(dir_name);
            end

            %make template for simulations
            template.vumps_opts = vumps_opts;
            template.handle = @make_PEPO_2D_A;

            template.dir_name = dir_name;
            template.time = dt;

            S_z = [1, 0; 0, -1];
            template.X = S_z; %observable

            
            template.pepo_opts = struct();

            saveboy(sprintf("%s/template.mat", template.name_prefix), 'template', template);
            
            first=1;
        end

        disp(template.name_prefix)

        template.vumps_opts.chi_max = chi;
        calc_ising_2d( 0.01, 0.01, nsammple, template, maxit,first);
    end
end

function calc_ising_2d(aim_dx, aim_dy, nsammple, template, maxit,first)

    for i = 1:maxit

        %determine next temperatures
        if first == 1
            
            x_min = template.x_bounds(1);
            x_max = template.x_bounds(2);
            
            x0 = (x_max - x_min) / (nsammple - 1) * (0:nsammple - 1) + x_min;
            first=0;
        else

            data = fetch_matfiles(template.name, struct);
            data = filter_ising_results(data, struct);

            x_arr = data.(  template.free_var );
            y_arr = data.m;

            dx_arr = diff(x_arr);
            dy_arr = diff(y_arr);

            if sum(abs(dx_arr) > aim_dx) == 0 && sum(abs(dy_arr) > aim_dy) == 0
                break;
            end

            ds = (dx_arr.^2 + ((aim_dx / aim_dy) .* dy_arr).^2).^(0.5);

            new_x = zeros(size(ds));

            for ii = 1:nsammple
                [~, idx] = max(ds ./ (new_x + 1));
                new_x(idx) = new_x(idx) + 1;
            end
            
            x0 = [];

            K = find(new_x ~= 0);
            for ii = 1:numel(K)
                ind = K(ii);
                nnx = new_x(ind);
                nx = x_arr(ind);
                dx = dx_arr(ind) / (nnx + 1);

                x0 = [x0, nx + dx * (1:nnx)];
            end

        end

        
        parfor j = 1:numel(x0)
            save_vars = [];
            save_vars.fname = sprintf('%2d:%2d', i, j);
            

            Ising2D_core(save_vars, template, x0(j));
        end
    end
end


