function vumpsObj = PEPO_vumps(pepo_matrix, vumps_opts, save_vars)

    if nargin < 3
        save_vars = [];
    end

    options = struct;

    if vumps_opts.disp == 1
        options.verbosity = Verbosity.Concise;
    else
        options.verbosity = Verbosity.None;
    end

    options.doPlot = true;
    options.doSave = false;
    options.tolInitial = 1e-4;
    options.tolConvergence = vumps_opts.tolfixed;
    options.dynamical = false;
    options.maxIterations = vumps_opts.vumps_maxit;
    options.schmidtCut = vumps_opts.tolfixed;

    options = VumpsOptions(options);

    %put into vumps format
    T = pepo_matrix;
    M = ncon({T}, {[1, 1, -1, -2, -3, -4]});


    if isfield(save_vars, 'G0') %legacy converter

        vumps_opts.chi_max = save_vars.A{1}.dims(1);
        
        o = TransferMpo( TensorNone(M) );
        
        initialMPS = UniformMps.FromCell( save_vars.A);
        

        initialConditions = struct('mps', initialMPS, 'environment', []);
        vumpsObj = Vumpser(o, initialConditions, options);
        
    elseif isfield(save_vars, 'vumpsObj')
        vumpsObj = save_vars.vumpsObj; %already correct, just finish calculation 
    else
        dims = [vumps_opts.chi_max, size(M, 1), vumps_opts.chi_max];
        
        
        cs = vumps_opts.cell_size;

        mps = UniformMps.empty;
        mpo = TensorNone.empty;

        for i = 1:cs(1)
            mps2 = TensorNone.empty;

            for j = 1:cs(2)
                mpo(i, j) = TensorNone(M);
                mps2(j) = TensorNone(rand(dims));
            end

            mps(i) = UniformMps(mps2);
        end

        o = TransferMpo(mpo);
        initialMPS = mps;
        
        initialConditions = struct('mps', initialMPS, 'environment', []);
        vumpsObj = Vumpser(o, initialConditions, options);
    end

    vumpsObj.DoVumps;

end
