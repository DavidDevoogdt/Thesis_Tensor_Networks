%function [mag, inv_corr_length, delta, ctr, err] = PEPO_get_expectation (obj, X, chimax, maxit, name, A, G0, T)


% 1x1: sampling_reproces('TIM_g=2.5_order_5_chi=20_sym=1_11_August_2021_16:24','V2');
% 2x2: sampling_reproces('TIM_g=2.5_order_5_chi=20_sym=1_11_August_2021_16:24','V2');
% 2x1: sampling_reproces('TIM_g=2.5_order_5_chi=15_sym=1_09_August_2021_10:42','V2');
% 1x2: sampling_reproces('TIM_g=2.5_order_5_chi=15_sym=1_12_August_2021_10:44_1x2','V2');
% 2x2: sampling_reproces('TIM_g=2.5_order_5_chi=15_sym=1_09_August_2021_10:00','V2');
% 2x2: sampling_reproces('TIM_g=2.5_order_5_chi=15_sym=1_12_August_2021_11:36_2x2','V2');  %(with extensions)

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

    skip_res = 0;%for debugging purposes

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

        
        rho2 = calculate_rho(vumpsObj, B, [m,o2],{1,2});
        mag2 = abs(trace(rho2 * X))

        rho2 = calculate_rho(vumpsObj, B, [o2],{1,1});
        mag2 = abs(trace(rho2 * X))
       
        %different
        rho2 = calculate_rho(vumpsObj, B, [o2,m],{1,1});
        mag2 = abs(trace(rho2 * X))
        
        
        tmag = 0;
        tmag2 = 0;
        
        for i = 1:vumpsObj.width
            for j=1:vumpsObj.depth
        
                rho = calculate_rho(vumpsObj, B, [o2,m;
                                                  m, m; ],{j,i});
                mag = abs(trace(rho * X))
                tmag = tmag+mag;
                
                
                rho2 = calculate_rho(vumpsObj, B, [o2],{j,i});
                mag2 = abs(trace(rho2 * X))
                tmag2 = tmag2+mag2;
            end
        end
        tmag = tmag/(vumpsObj.width*vumpsObj.width);
        tmag2 = tmag2/(vumpsObj.width*vumpsObj.width);


        sv2 = MpsSchmidtValues( vumpsObj.mps(1).C(1)  ).^2;
        S = -sum(sv2 .* log(sv2));


            %calculate epsilon_i
            if p.Results.doEpsi || p.Results.doVumps
                opts = [];
                opts.krylovdim = 100; opts.tol = 1e-14;
                opts.level = 1;
        
                %[~, f] = TensorEigs(@(x) get_Bc(x), B, 8, 'lm', opts);
                [~, f] = TensorEigs(@(x) get_Ac(x), Ac, 8, 'lm', opts);
                f2 = f(2:end) ./ f(1);
        
                eps_i = -log(f2);
                inv_corr_length = eps_i(1);
        
                delta = eps_i(4) - eps_i(2);
        
                results.marek = delta;
                results.eps_i = eps_i;
        
                results.inv_corr_length = inv_corr_length;
        
            end

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
