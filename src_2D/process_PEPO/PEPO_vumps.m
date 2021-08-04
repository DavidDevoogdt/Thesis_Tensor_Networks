function [A, G1, lambda1, ctr, err] = PEPO_vumps(pepo_matrix, vumps_opts, save_vars)

    if nargin < 3
        save_vars = [];
    end

    %     opts.charges = 'regular';
    %     opts.dynamical = 'off';
    %     opts.dyncharges = 0;
    %     opts.schmidtcut = 1e-10;
    %     opts.chimax = 350;
    %     opts.disp = vumps_opts.disp;
    %     opts.tolmax = 1e-4;
    %     opts.tolfactor = 1e4;
    %     opts.minit = 10;
    %     opts.dyniter = 5;
    %     opts.truncate = 0;
    %     opts.method = 'vumps';
    %
    %     opts.plot = 'on';
    %     opts.maxit = vumps_opts.vumps_maxit;
    %     opts.tolfixed = vumps_opts.tolfixed;

    options = struct;
    options.verbosity = Verbosity.Concise;
    options.doPlot = true;
    options.doSave = false;
    options.tolInitial = 1e-4;
    options.tolConvergence = 1e-5;
    options.dynamical = false;

    options = VumpsOptions(options);

    %put into vumps format
    T = pepo_matrix;
    M = ncon({T}, {[1, 1, -1, -2, -3, -4]});

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

    if isfield(save_vars, 'G0')
        error('todo');
        [A, G1, lambda1, ctr, err] = Vumps(o, save_vars.A, save_vars.G0, opts);
    else
        initialConditions = struct('mps', initialMPS, 'environment', []);
        vumpsObj = Vumpser(o, initialConditions, options);
        vumpsObj.DoVumps;

        %[A, G1, lambda1, ctr, err] = Vumps(o, vumps_opts.chi_max, [], opts);
    end

end

function a = get_mps(dims)
    a = TensorNone(rand(dims));
end
