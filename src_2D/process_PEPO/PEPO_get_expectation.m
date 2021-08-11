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
        vumpsObj = PEPO_vumps(T, vumps_opts, save_vars);

        results.ctr = numel(vumpsObj.history.errors);
        results.err = vumpsObj.error;
    else
        vumpsObj = save_vars.vumpsObj;
    end

    skip_res = 0;

    if skip_res == 0

        %construct central tensor
        m = get_tensor(T);
        %mx = get_tensor(T, X);

        o2 = TensorNone(ncon({T}, {[-5, -6, -1, -2, -3, -4]}));

        doSym = 1;
        if doSym
            B = vumpsObj.mps; %mps for rotated system is identical
        else
            error('check')

        end

        
        for i=1:vumpsObj.width
            for j=1:vumpsObj.depth
        
                [rho, f, lambda] = calculate_rho(vumpsObj, B, [o2],{j,i});
                mag = abs(trace(rho * X))
                
                
                [rho2,f2,l]= calculate_rho(vumpsObj, B, [m,o2;
                                                         m,m ],{j,i-1});
                mag = abs(trace(rho2 * X))
                
                [rho2,f2,l]= calculate_rho(vumpsObj, B, [m,o2;],{j,i-1});
                mag = abs(trace(rho2 * X))
                
            end
        end
      

        
        [rho2,f2,ll]= calculate_rho(vumpsObj, B, [o2,m]);

        %[rho3,f3] = calculate_rho(A,B,G0,[o2,m,o2]);

        %rho3 = calculate_rho(A,B,G0,F,[m; m],X); % too expensive
        %rho4 = calculate_rho(A,B,G0,F,[m, o2;
        %                               m, m],X);

        mag = abs(trace(rho * X));

        %get entanglement entropy

        sv2 = MpsSchmidtValues(Ca).^2;
        S = -sum(sv2 .* log(sv2));

        %S =  -trace( rho * logm(rho)) ; %doesn't work, not enough
        %entanglement possible
        %
        %         else
        %             %error('todo')
        %             mag = 0;
        %             S = 0;
        %
        %         end

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
    else
        mag = 0;
        S = 0;
        chi = 0;
    end

    save_vars.vumpsObj = vumpsObj;
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

    m = TensorNone(M);
    
end
