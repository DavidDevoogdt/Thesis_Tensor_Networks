function succes_string = sampling_core(save_vars, template, x, results, get_exp_opts)
    %core routine: calculates pepo according to template and temperature/transversal field x,
    %calculates environment with PEPO and saves to disk in 2 files: results (lightweight) and save_vars (Vumps environment,...)
    %depending on save_vars, just recalculates the expectation value with previous environment and PEPO

    if nargin < 5
        get_exp_opts = struct();
    end

    if nargin < 4
        results = struct;
    end

    switch template.fixed_var
        case 'g'
            save_vars.model_params = models(template.model_name, struct('g', template.fixed_val));
            beta = 1 / x;
        case 'T'
            beta = 1 / template.fixed_val;
            save_vars.model_params = models('t_ising', struct('g', x));
    end

    if ~isfield(save_vars, 'PEPO_matrix')

        popts = template.pepo_opts;
        popts.beta = beta;

        pepo = PEPO(save_vars.model_params, popts, template.handle);

        if pepo.error_code ~= 0
            fprintf("%s %3d:%s %s:%.4e: construction failed\n", datestr(now, 'HH:MM:SS'), template.vumps_opts.chi_max, strrep(save_vars.fname, '_', ':'), template.free_var, x)
            return;
        end

        save_vars.PEPO_matrix = pepo.PEPO_matrix;

        save_vars.virtual_level_sizes_horiz = pepo.virtual_level_sizes_horiz;
        save_vars.virtual_level_sizes_vert = pepo.virtual_level_sizes_vert;

    end

    [results, save_vars] = PEPO_get_expectation (template.X, save_vars, template.vumps_opts, results, get_exp_opts);

    results.(template.fixed_var) = template.fixed_val;
    results.(template.free_var) = x;

    save_vars.results_name = sprintf("%s/results_%s.mat", template.name_prefix, save_vars.fname);
    save_vars.save_vars_name = strrep(save_vars.results_name, '/results_', '/save_vars_');

    saveboy(save_vars.results_name, 'results', results);
    saveboy(save_vars.save_vars_name, 'save_vars', save_vars);

    succes_string = sprintf("%s %3d:%s %s:%.4e <X>:%.4e xi:%.4e marek gap:%.4f ctr:%3d err:%.4e\n", datestr(now, 'HH:MM:SS'), results.chi, strrep(save_vars.fname, '_', ':'), template.free_var, x, results.m, 1 / results.inv_corr_length, results.marek, results.ctr, results.err);

end
