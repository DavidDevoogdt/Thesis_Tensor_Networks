%function [mag, inv_corr_length, delta, ctr, err] = PEPO_get_expectation (obj, X, chimax, maxit, name, A, G0, T)

function [results, save_vars] = PEPO_get_expectation (X, save_vars, vumps_opts,results,opts)

    if nargin<5
       opts=[]; 
    end
    
    p = inputParser;
    addParameter(p, 'doVumps', ~isfield(save_vars, 'A') || ~isfield(save_vars, 'G0') || opts.doVumps)
    addParameter(p, 'doEpsi', 1)

    parse(p, opts)
    
    assert(isfield(save_vars, 'PEPO_matrix'));
    T = save_vars.PEPO_matrix;

    if p.Results.doVumps
        [A, G0, ~, ctr, err] = PEPO_vumps(T, vumps_opts,save_vars);
        
        results.ctr = ctr;
        results.err = err;
    else
        A = save_vars.A;
        G0 = save_vars.G0;
    end

    %construct central tensor
    M = ncon({T}, {[1, 1, -1, -2, -3, -4]});
    m.legs = 4;
    m.group = 'none';
    m.dims = size(M);
    m.var = M;
    
    %same but open
    O2= ncon( {T},{[-5,-6,-1,-2,-3,-4]}  ) ;
    o2.legs = 6;
    o2.group = 'none';
    o2.dims = size(O2);
    o2.var = O2;
    
    % calculate eigentensor under
    GL = G0{1}; GR = G0{2}; Ac = A{4};

    %get Ac equivalent from under
    opts = [];
    opts.krylovdim = 100; opts.tol = 1e-14;
    opts.level = 1;

    function x = transfer_down(x)
        x = TensorContract({x, GL, m, GR}, ...
            {[1, 2, 5], [1, 3, -1], [-2, 4, 2, 3], [-3, 4, 5]});
    end

    if p.Results.doVumps || ~isfield(save_vars, 'B')
        [B, ~, ~] = TensorEigs(@(x) transfer_down(x), TensorConj(Ac), 1, 'lm', opts);
    else
        B = save_vars.B;
    end

    %calculate <X>
    
    [rho,~] = TensorContract({B, GL, Ac, o2, GR}, ...
        {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3,-1,-2], [8, 7, 6]});
   
    mag = abs(  trace(rho.var*X) /  trace(rho.var) );

    %get entanglement entropy
    
    sv = MpsSchmidtValues(A).^2;
    
    S = -sum( sv.*log(sv)  );
    
    %calculate epsilon_i
    if p.Results.doEpsi || p.Results.doVumps 
        opts = [];
        opts.krylovdim = 100; opts.tol = 1e-14;
        opts.level = 1;

        [~, f] = TensorEigs(@(x) transfer_down(x), B, 8, 'lm', opts);
        f2 = f(2:end) ./ f(1);

        eps_i = -log(f2);
        inv_corr_length = eps_i(1);

        delta = eps_i(4) - eps_i(2);
        
        results.marek = delta;
        results.eps_i = eps_i;
          
        results.inv_corr_length = inv_corr_length;

    end

    %assign output
    save_vars.A = A;
    save_vars.G0 = G0;
    save_vars.B = B;
    save_vars.PEPO_matrix = T;

    results.m = mag;
    
    
    results.S = S;
    results.chi = vumps_opts.chi_max;
    
    
    results.ftime = now;

end
