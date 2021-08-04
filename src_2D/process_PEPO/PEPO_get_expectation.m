%function [mag, inv_corr_length, delta, ctr, err] = PEPO_get_expectation (obj, X, chimax, maxit, name, A, G0, T)

function [results, save_vars] = PEPO_get_expectation (X, save_vars, vumps_opts, results, opts)

    if nargin < 5
        opts = [];
    end

    p = inputParser;
    addParameter(p, 'doVumps', ~isfield(save_vars, 'A') ||~isfield(save_vars, 'G0') || opts.doVumps)
    addParameter(p, 'doEpsi', 1)

    parse(p, opts)

    assert(isfield(save_vars, 'PEPO_matrix'));
    T = save_vars.PEPO_matrix;

    if p.Results.doVumps
        [A, G0, ~, ctr, err] = PEPO_vumps(T, vumps_opts, save_vars);

        results.ctr = ctr;
        results.err = err;
    else
        A = save_vars.A;
        G0 = save_vars.G0;
    end

    %construct central tensor
    m = get_tensor(T);
    mx = get_tensor(T, X);

    %same but open
    O2 = ncon({T}, {[-5, -6, -1, -2, -3, -4]});
    o2.legs = 6;
    o2.group = 'none';
    o2.dims = size(O2);
    o2.var = O2;

    if vumps_opts.cell_size == 1

        [GL, Gr] = G0{:};
        [Al, Ar, Ca, Ac] = A{:};
        chi = Ar.dims(1);

        doSym = 1;
        if doSym

            Bl = Ar; Bl.var = ncon({Bl.var}, {[-3, -2, -1]});
            Br = Al; Br.var = ncon({Br.var}, {[-3, -2, -1]});
            Bc = Ac; Bc.var = ncon({Bc.var}, {[-3, -2, -1]});
            Cb = Ca; Cb.var = ncon({Cb.var}, {[-2, -1]});
            B = {Bl, Br, Cb, Bc};

        else
            error('check')
            [v, ~, ~] = TensorEigs(@(x) get_Bc(x, GL, GR), TensorConj(Ac), 1, 'lm', opts);
        end

        %[rho,f]  = calculate_rho(A,B,G0,[o2]);

        %[rho2,f2]= calculate_rho(A,B,G0,[m,o2]);

        %[rho3,f3] = calculate_rho(A,B,G0,[o2,m,o2]);

        %rho3 = calculate_rho(A,B,G0,F,[m; m],X); % too expensive
        %rho4 = calculate_rho(A,B,G0,F,[m, o2;
        %                               m, m],X);

        mag = abs(trace(rho.var * X));

        %get entanglement entropy

        sv2 = MpsSchmidtValues(Ca).^2;
        S = -sum(sv2 .* log(sv2));

        %S =  -trace( rho * logm(rho)) ; %doesn't work, not enough
        %entanglement possible

    else
        %error('todo')
        mag = 0;
        S = 0;

    end

    %     function x = get_Ac(x)
    %         x = TensorContract({GL, x, GR, m}, ...
    %             {[-1, 2, 1], [1, 3, 4], [4, 5, -3], [3, 5, -2, 2]});
    %     end
    %
    %     function x = get_C(x)
    %         x = TensorContract({GL, x, GR}, {[-1, 3, 1], [1, 2], [2, 3, -2]});
    %     end
    %
    %     %
    %     function x = get_Bc(x, GL, GR)
    %         x = TensorContract({x, GL, m, GR}, ...
    %             {[1, 2, 5], [1, 3, -1], [-2, 4, 2, 3], [-3, 4, 5]});
    %     end
    %
    %     function x = get_C_down(x)
    %         x = TensorContract({x, GL, GR}, ...
    %             {[1, 3], [1, 2, -1], [-2, 2, 3]});
    %     end
    %
    %        function a = overlap(G, A)
    %         a = TensorContract({G{1}, A{3}, G{2}, TensorConj(A{3})}, ...
    %             {[1, 4:A{1}.legs + 1, 2], [2, 3], [3, 4:A{1}.legs + 2], [1, A{1}.legs + 2]});
    %     end
    %
    %
    %     function x = MVumpsGL_ApplyFunction1(x, d, AL, cAL, O)
    %         dd = mod(d, depth) + 1;
    %         [~, width] = size(AL);
    %         for ww = 1:width
    %             x = TensorContract({cAL{dd, ww}, x, AL{d, ww}, O}, ...
    %                 {[1, 2, -1], [1, 3, 4], [4, 5, -3], [5, -2, 2, 3]});
    %         end
    %     end
    %
    %
    %     %calculate epsilon_i
    %     if p.Results.doEpsi || p.Results.doVumps
    %         opts = [];
    %         opts.krylovdim = 100; opts.tol = 1e-14;
    %         opts.level = 1;
    %
    %         %[~, f] = TensorEigs(@(x) get_Bc(x), B, 8, 'lm', opts);
    %         [~, f] = TensorEigs(@(x) get_Ac(x), Ac, 8, 'lm', opts);
    %         f2 = f(2:end) ./ f(1);
    %
    %         eps_i = -log(f2);
    %         inv_corr_length = eps_i(1);
    %
    %         delta = eps_i(4) - eps_i(2);
    %
    %         results.marek = delta;
    %         results.eps_i = eps_i;
    %
    %         results.inv_corr_length = inv_corr_length;
    %
    %     end

    %assign output
    save_vars.A = A;
    save_vars.G0 = G0;
    save_vars.PEPO_matrix = T;

    results.m = mag;

    results.S = S;
    results.chi = chi;

    results.ftime = now;

end

function m = get_tensor(T, X)

    if nargin < 2
        X = eye(2);
    end

    M = ncon({T, X}, {[1, 2, -1, -2, -3, -4], [2, 1]});

    m.legs = 4;
    m.group = 'none';
    m.dims = size(M);
    m.var = M;
end
