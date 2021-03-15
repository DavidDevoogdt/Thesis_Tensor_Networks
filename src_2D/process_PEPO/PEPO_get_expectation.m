function [mag, inv_corr_length, delta, ctr, err] = PEPO_get_expectation (obj, X, chimax, maxit, name, A, G0, ctr, err, T)
    if nargin <= 5%calculate stuff
        [A, G0, ~, ctr, err] = PEPO_vumps(obj, chimax, maxit, name);
        T = obj.PEPO_matrix;
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
    opts.krylovdim = 100; opts.tol = 1e-14;

    function x = transfer_down(x)
        x = TensorContract({x, GL, m, GR}, ...
            {[1, 2, 5], [1, 3, -1], [-2, 4, 2, 3], [-3, 4, 5]});
    end

    [B, ~, ~] = TensorEigs(@(x) transfer_down(x), TensorConj(Ac), 1, 'lm', opts);

    %calculate <X>
    [x, ~] = TensorContract({B, GL, Ac, m, GR}, ...
        {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3], [8, 7, 6]});

    [y, ~] = TensorContract({B, GL, Ac, o, GR}, ...
        {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3], [8, 7, 6]});

    mag = y / x;

    %calculate transfer matrix eigenvalues and extract xi and delta
    function x = transfer_up(x)
        [x, ~] = TensorContract({GL, x, m, GR}, ...
            {[-1, 3, 4], [4, 5, 8], [5, 7, -2, 3], [8, 7, -3]});
    end

    [~, f] = TensorEigs(@(x) transfer_up(x), Ac, 10, 'lm', opts);

    f2 = f(2:end) ./ f(1);

    eps_i = -log(f2);
    inv_corr_length = eps_i(1);

    delta = eps_i(2) - eps_i(1);

end
