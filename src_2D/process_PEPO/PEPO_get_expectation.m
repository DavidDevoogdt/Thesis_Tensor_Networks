%function [mag, inv_corr_length, delta, ctr, err] = PEPO_get_expectation (obj, X, chimax, maxit, name, A, G0, T)


% 1x1: sampling_reproces('TIM_g=2.5_order_5_chi=20_sym=1_11_August_2021_16:24','V2');
% 2x2: sampling_reproces('TIM_g=2.5_order_5_chi=20_sym=1_11_August_2021_16:24','V2');
% 2x1: sampling_reproces('TIM_g=2.5_order_5_chi=15_sym=1_09_August_2021_10:42','V2');
% 1x2: sampling_reproces('TIM_g=2.5_order_5_chi=15_sym=1_12_August_2021_10:44_1x2','V2');
% 2x2: sampling_reproces('TIM_g=2.5_order_5_chi=15_sym=1_09_August_2021_10:00','V2');
% 2x2: sampling_reproces('TIM_g=2.5_order_5_chi=10_sym=1_12_August_2021_13:37_2x2','V2');  %(with extensions)


function [results, save_vars] = PEPO_get_expectation (X, save_vars, vumps_opts, results, opts)

    if nargin < 5
        opts = [];
    end

  
    p = inputParser;
    addParameter(p, 'doVumps', ~isfield(save_vars, 'vumpsObj')  )
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
            error('todo')

        end

        rho = calculate_rho(vumpsObj, B, [o2]);
        mag = abs(trace(rho * X));

        sv2 = MpsSchmidtValues( vumpsObj.mps(1).C(1)  ).^2;
        S = -sum(sv2 .* log(sv2));


        %calculate epsilon_i
        if p.Results.doEpsi || p.Results.doVumps
           
            f = TransferEigs( vumpsObj.mps(1), vumpsObj.mps(1),[], 7 );
            
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
    end

    if isfield( save_vars,'A')
        save_vars = rmfield(save_vars, 'A');
        save_vars = rmfield(save_vars, 'G0');
    end
    
    save_vars.vumpsObj = vumpsObj;
    save_vars.PEPO_matrix = T;

    results.m = mag;
    results.S = S;
    results.chi = vumpsObj.mps(1).bondDimension;
    results.ftime = now;

end

function m = get_tensor(T, X)

    if nargin < 2
        X = eye(2);
    end

    M = ncon({T, X}, {[1, 2, -1, -2, -3, -4], [2, 1]});

    m = TensorNone(M);
    
end
