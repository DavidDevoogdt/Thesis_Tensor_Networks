function Ising2D_core(save_vars, template, x, results,get_exp_opts)

    if nargin < 5
        get_exp_opts = struct();
    end

    if nargin < 4
        results = [];
    end

    switch template.fixed_var
        case 'g'
            save_vars.model_params = models( template.model_name  , struct('g', template.fixed_val ));
            beta = 1/x;
        case 'T'
            beta = 1/template.fixed_val;
            save_vars.model_params = models('t_ising', struct('g', x));
    end
    
    
    


    if ~isfield(save_vars, 'PEPO_matrix')

        pepo = PEPO(save_vars.model_params.d, -beta * save_vars.model_params.H_1_tensor, -beta * save_vars.model_params.H_2_tensor, 5, template.handle, template.pepo_opts);
        %save_vars = [];
        save_vars.PEPO_matrix = pepo.PEPO_matrix;
        
        save_vars.virtual_level_sizes_horiz = pepo.virtual_level_sizes_horiz;
        save_vars.virtual_level_sizes_vert = pepo.virtual_level_sizes_vert;

    end
    
    [results, save_vars] = PEPO_get_expectation (template.X, save_vars, template.vumps_opts, results,get_exp_opts);
    
    

    results.( template.fixed_var ) =  template.fixed_val;
    results.( template.free_var  ) = x;
     

    save_vars.results_name = sprintf("%s/results_%s.mat", template.name_prefix, save_vars.fname);
    save_vars.save_vars_name = strrep(save_vars.results_name, '/results_', '/save_vars_');

    saveboy(save_vars.results_name, 'results', results);
    saveboy(save_vars.save_vars_name, 'save_vars', save_vars);

    fprintf("%s %3d:%s %s:%.4e <X>:%.4e xi:%.4e marek gap:%.4f ctr:%3d err:%.4e\n", datestr(now, 'HH:MM:SS'), template.vumps_opts.chi_max, strrep(save_vars.fname, '_', ':') ,  template.free_var , x   , results.m, 1 / results.inv_corr_length, results.marek, results.ctr, results.err);
end
