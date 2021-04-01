function Ising2D_core(save_vars, template, T, results)

    if nargin < 4
        results = [];
    end

    if ~isfield(save_vars, 'PEPO_matrix')

        beta = 1 / T;

        pepo = PEPO(template.model_params.d, -beta * template.model_params.H_1_tensor, -beta * template.model_params.H_2_tensor, 5, template.handle, template.pepo_opts);
        %save_vars = [];
        save_vars.PEPO_matrix = pepo.PEPO_matrix;

    end

    [results, save_vars] = PEPO_get_expectation (template.X, save_vars, template.vumps_opts, results);
    results.T = T;

    save_vars.results_name = sprintf("%s/results_%s.mat", template.name_prefix, save_vars.fname);
    save_vars.save_vars_name = strrep(save_vars.results_name, '/results_', '/save_vars_');

    saveboy(save_vars.results_name, 'results', results);
    saveboy(save_vars.save_vars_name, 'save_vars', save_vars);

    fprintf("%s %3d:%s T:%.4e mag:%.4e xi:%.4e marek gap:%.4f ctr:%3d err:%.4e\n", datestr(now, 'HH:MM:SS'), template.vumps_opts.chi_max, strrep(save_vars.fname, '_', ':'), results.T, results.m, 1 / results.inv_corr_length, results.marek, results.ctr, results.err);

end

function m = m_onsager(T, J)
    T_c = 2 * J / (log(1 + sqrt(2)));
    m = T;
    mask = T < T_c;
    m(mask) = (1 - sinh((2 * J) ./ T(mask)).^(-4)).^(1/8);
    m(~mask) = 0;
end
