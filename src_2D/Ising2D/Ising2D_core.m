function Ising2D_core(name_prefix, T0,template,ind)

    beta = 1 ./ T0;

    %dir for temp results
    dir_name = sprintf("%s/", name_prefix);
    if ~exist(dir_name, 'dir')
        mkdir(dir_name);
    end

    %hamiltonian setup
    S_x = [0, 1; 1, 0];
    S_y = [0, -1i; 1i, 0];
    S_z = [1, 0; 0, -1];
    I_tensor = eye(2);

    handle = @make_PEPO_2D_A;

    d = 2;

    J=template.J;
    g=template.g;
    
    template.X=S_z;
    
    H_1_tensor = -J * g * S_x;
    H_2_tensor = -J * (reshape(ncon({S_z, S_z}, {[-1, -3], [-2, -4]}), [d, d, d, d]));

    opts = [];
    opts.testing = 0;
    opts.visualise = 0;
    opts.double = 0;

    parfor iii = 1:numel(T0)

            pepo = PEPO(d, -beta(iii) * H_1_tensor, -beta(iii) * H_2_tensor, 5, handle, opts);
         
            arr = pepo.PEPO_matrix;
            
            in_vars = [];
            in_vars.PEPO_matrix = arr;
            
            vumps_opts=[];
            vumps_opts.chi_max = template.chi;
            vumps_opts.maxit = template.vumps_maxit;
            
            [results,save_vars] = PEPO_get_expectation (S_z, in_vars, vumps_opts );
            results.T= T0(iii);
            
            saveboy(sprintf("%s/results_%d_%d.mat", name_prefix,ind, iii), 'results','template',results,template);
            saveboy(sprintf("%s/save_vars_%d_%d.mat", name_prefix,ind, iii), 'save_vars',save_vars);
            
           fprintf("%s %3d:%2d:%2d T:%.4e mag:%.4e xi:%.4e marek gap:%.4f ctr:%3d err:%.4e\n", datestr(now, 'HH:MM:SS'), template.chi, ind, iii, results.T, results.m, 1/ results.inv_corr_length, results.marek, results.ctr, results.err);

    end

end

function m = m_onsager(T, J)
    T_c = 2 * J / (log(1 + sqrt(2)));
    m = T;
    mask = T < T_c;
    m(mask) = (1 - sinh((2 * J) ./ T(mask)).^(-4)).^(1/8);
    m(~mask) = 0;
end
