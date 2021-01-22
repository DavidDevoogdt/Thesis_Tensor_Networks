function [mag, corr_length, delta] = PEPO_get_expectation (obj, X, chimax)
    [A, B, G, lambda] = vumps(obj, chimax);

    T = obj.PEPO_matrix;

    M = ncon({T}, {[1, 1, -1, -2, -3, -4]});

    m.legs = 4;
    m.group = 'none';
    m.dims = size(M);
    m.var = M;

    O = ncon ({T, X}, {[1, 2, -1, -2, -3, -4], [1, 2]});
    o.legs = 4;
    o.group = 'none';
    o.dims = size(O);
    o.var = O;

    %transfereigs edited !

    GL = G{1}; GR = G{2}; Ac = A{4};

    %should be lambda?
    [x, ~] = TensorContract({B, GL, Ac, m, GR}, ...
        {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3], [8, 7, 6]});

    [y, ~] = TensorContract({B, GL, Ac, o, GR}, ...
        {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3], [8, 7, 6]});

    mag = y / x;

    function x = transfer_up(x)
        [x, ~] = TensorContract({GL, x, m, GR}, ...
            {[-1, 3, 4], [4, 5, 8], [5, 7, -2, 3], [8, 7, -3]});
    end

    opts.krylovdim = 100; opts.tol = 1e-14;

    [rho, f] = TensorEigs(@(x) transfer_up(x), A{4}, 5, 'lm', opts);

    eps_i = -log(abs(f));
    corr_length = eps_i(1 + 1) - eps_i(1);

    delta = eps_i(4 + 1) - eps_i(2 + 1); %same as https://arxiv.org/pdf/1907.08603.pdf
end
