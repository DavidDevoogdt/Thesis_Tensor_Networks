
%test_2D(struct('do_loops',0,'loop_extension',0,'inv_eps',1e-10))

function test_2D( opts,map  )
    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');

    pathparts = [pathparts(1:end - 3), 'test_2D_files'];
    fold2 = strjoin(pathparts, '/');
    time = now;

    filename = sprintf("%s/2D_%s.mat", fold2, datestr(time, 'mm-dd-yy_HH-MM-SS'));
    fprintf("%s \n", filename)

    model_opts.g = 2.5;
    model = 't_ising';
    simul = models(model, model_opts);

  
    handle = @make_PEPO_2D_A;

    beta_arr = 10.^(-3 : 0.1 : 2);
    beta_len = size(beta_arr, 2);
    err_arr = zeros(beta_len, 1);

    
    if map == 1
    
        num_map = [
            11,0,0,0,0,0;
            10,0,0,0,0,0;
            9,0,0,0,0,0;
            8,0,0,0,0,0;
            7,12,13,0,0,0;
            1,2,3,4,5,6;];
        map_opts = struct("numbered", true, "h_cyclic", 1, "v_cyclic", 1);
        density_site = 1;
        
        map = create_map(num_map,map_opts);
        
    else
        num_map = [
            0, 0, 13, 14,  0, 0;
            7, 8,  9, 10, 11, 12;
            1, 2,  3,  4,  5, 6;];
        map_opts = struct("numbered", true, "h_cyclic", 0, "v_cyclic", 0);
        density_site = 9;
        
        map = create_map(num_map,map_opts);
    end
    
    %parfor i = 1:beta_len
    parfor i = 1:beta_len
        opts2 = opts;
        opts2.beta = beta_arr(i);
        try
            [pepo,err_code] = PEPO(simul, opts2, handle);
            if err_code ==0
                err = calculate_error(pepo, map, [], 1, density_site);
            else
               err = NaN; 
            end
        catch
            err = NaN;
        end

        err_arr(i) = abs(err); 
        fprintf(" beta %.4e cycl err %.4e \n", beta_arr(i), err);
    end
    
    saveboy(filename, 'time', 'simul', 'beta_arr', 'err_arr', 'num_map', 'map_opts', 'density_site','opts',time, simul, beta_arr, err_arr, num_map, map_opts, density_site,opts)
end
