function [A, G1, lambda1, ctr, err] = PEPO_vumps(pepo_matrix, vumps_opts, save_vars)

    if nargin < 3
        save_vars = [];
    end

    opts.charges = 'regular';
    opts.dynamical = 'off';
    opts.dyncharges = 0;
    opts.schmidtcut = 1e-10;
    opts.chimax = 350;
    opts.disp = vumps_opts.disp;
    opts.tolmax = 1e-4;
    opts.tolfactor = 1e4;
    opts.minit = 10;
    opts.dyniter = 5;
    opts.truncate = 0;
    opts.method = 'vumps';

    opts.plot = 'on';
    opts.maxit = vumps_opts.vumps_maxit;
    opts.tolfixed = vumps_opts.tolfixed;

    %put into vumps format
    T = pepo_matrix;
    M = ncon({T}, {[1, 1, -1, -2, -3, -4]});

    o.legs = 4;
    o.group = 'none';
    o.dims = size(M);
    o.var = M;

    if vumps_opts.cell_size == 1
        O.type = 'mpo';
        O.mpo = o;
    else
        O.type = 'multi_mpo'; O.mpo = {o, o; o, o};
    end

    if isfield(save_vars, 'G0')
        [A, G1, lambda1, ctr, err] = Vumps(O, save_vars.A, save_vars.G0, opts);
    else
        [A, G1, lambda1, ctr, err] = Vumps(O, vumps_opts.chi_max, [], opts);
    end

end
