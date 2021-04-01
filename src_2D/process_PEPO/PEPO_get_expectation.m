%function [mag, inv_corr_length, delta, ctr, err] = PEPO_get_expectation (obj, X, chimax, maxit, name, A, G0, T)

function [results, save_vars] = PEPO_get_expectation (X, save_vars, vumps_opts,results,opts)

    if ~isfield(opts, 'doVumps')
        opts.doVumps = 1;
    end

    assert(isfield(save_vars, 'PEPO_matrix'));
    T = save_vars.PEPO_matrix;

    if ~isfield(save_vars, 'A') || ~isfield(save_vars, 'G0') || opts.doVumps
        [A, G0, ~, ctr, err] = PEPO_vumps(T, vumps_opts);
    else
        A = save_vars.A;
        G0 = save_vars.G0;
        ctr = results.ctr;
        err = results.err;
    end

    %construct central tensor
    M = ncon({T}, {[1, 1, -1, -2, -3, -4]});
    m.legs = 4;
    m.group = 'none';
    m.dims = size(M);
    m.var = M;

    %same but with operator
    O = ncon ({T, X}, {[1, 2, -1, -2, -3, -4], [1, 2]});
    o.legs = 4;
    o.group = 'none';
    o.dims = size(O);
    o.var = O;

    % calculate eigentensor under
    GL = G0{1}; GR = G0{2}; Ac = A{4};

    %get Ac equivalent from under
    opts = [];
    opts.krylovdim = 100; opts.tol = 1e-14;
    opts.level = 1;
    %opts.disp = 'iter-detailed';

    function x = transfer_down(x)
        x = TensorContract({x, GL, m, GR}, ...
            {[1, 2, 5], [1, 3, -1], [-2, 4, 2, 3], [-3, 4, 5]});
    end

    if ~isfield(save_vars, 'B')
        [B, ~, ~] = TensorEigs(@(x) transfer_down(x), TensorConj(Ac), 1, 'lm', opts);
    else
        B = save_vars.B;
    end

    %calculate <X>
    [x, ~] = TensorContract({B, GL, Ac, m, GR}, ...
        {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3], [8, 7, 6]});

    [y, ~] = TensorContract({B, GL, Ac, o, GR}, ...
        {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3], [8, 7, 6]});

    mag = abs(y / x);

    %calculate transfer matrix eigenvalues and extract xi and delta
    %     function x = transfer_up(x)
    %         [x, ~] = TensorContract({GL, x, m, GR}, ...
    %             {[-1, 3, 4], [4, 5, 8], [5, 7, -2, 3], [8, 7, -3]});
    %     end

    opts = [];
    opts.krylovdim = 100; opts.tol = 1e-14;
    opts.level = 1;

    [~, f] = TensorEigs(@(x) transfer_down(x), B, 8, 'lm', opts);

    %[~, f] = TensorEigs(@(x) transfer_down(x), B, 5, 'lm', opts);

    f2 = f(2:end) ./ f(1);

    eps_i = -log(f2);
    inv_corr_length = eps_i(1);

    delta = eps_i(4) - eps_i(2);

    s = MpsEntSpectrum(A);
    s2 = s .* conj(s);
    S = -sum(s2 .* log(s2));

    %assign output
    save_vars.A = A;
    save_vars.G0 = G0;
    save_vars.B = B;
    save_vars.PEPO_matrix = T;

    results.m = mag;
    results.marek = delta;
    results.eps_i = eps_i;
    results.ctr = ctr;
    results.err = err;
    results.inv_corr_length = inv_corr_length;
    results.S = S;
    
    results.ftime = now;

end
