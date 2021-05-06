function test_2D(model)
    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');

    pathparts = [pathparts(1:end - 3), 'test_2D_files'];
    fold2 = strjoin(pathparts, '/');

    time = now;

    filename = sprintf("%s/2D_%s.mat", fold2, datestr(time, 'mm-dd-yy_HH-MM-SS'));

    fprintf("%s \n", filename)

    model_opts.g = 2.5;
    %model = 't_ising';
    model = 'Heisenberg_2D';
    %model = 'Heisenberg_2D_X';

    simul = models(model, model_opts);

    opts.testing = 1;
    opts.visualise = 0;
    opts.double = 0;
    opts.inv_eps = 1e-12;
    opts.err_tol = 1e-13;

    pepo_order = 6;
    handle = @make_PEPO_2D_A;

    beta_arr = 10.^(log10(1):0.5:2);

    beta_len = size(beta_arr, 2);
    err_arr = zeros(beta_len, 1);

    num_map = [
            0, 5, 6, 0;
            1, 4, 7, 10;
            2, 3, 8, 9; ];
    map_opts = struct("numbered", true, "h_cyclic", 1, "v_cyclic", 0);
    density_site = 6;

    for i = 1:beta_len
        beta = beta_arr(i);

        pepo = PEPO(simul.d, -beta * simul.H_1_tensor, ...
            -beta * simul.H_2_tensor, ...
            pepo_order, handle, opts);

        err = calculate_error(pepo, num_map, map_opts, 1, density_site);

        fprintf(" beta %.4e cycl err %.4e \n", beta, err);
        err_arr(i) = abs(err);

        save(filename, 'time', 'simul', 'beta_arr', 'err_arr', 'num_map', 'map_opts', 'density_site')
    end
end
