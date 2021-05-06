function [A, G1, lambda1, ctr, err] = PEPO_vumps(pepo_matrix, vumps_opts, save_vars)

    if nargin < 3
        save_vars = [];
    end

    %     p = inputParser;
    %
    %     addRequired(p, 'chi_max')
    %     addParameter(p, 'name',[])
    %     addRequired(p, 'maxit' )
    %
    %     parse(p, vumps_opts)

    %todo check these params
    opts.charges = 'regular';
    opts.dynamical = 'off';
    opts.dyncharges = 0;
    opts.schmidtcut = 1e-10;
    opts.chimax = 350;
    %opts.disp='iter';
    opts.disp = 'none';
    %opts.disp = 'conv';
    opts.tolmax = 1e-4; %1e-4
    opts.tolfactor = 1e4;
    opts.minit = 1;
    opts.dyniter = 5;
    opts.truncate = 0;
    opts.method = 'vumps';

    %opts.save = 1;
    %opts.name = p.Results.name;

    %opts.method = 'qr';

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

    O.type = 'mpo';
    O.mpo = o;

    if isfield(save_vars, 'G0')
        [A, G1, lambda1, ctr, err] = Vumps(O, save_vars.A, save_vars.G0, opts);
    else
        [A, G1, lambda1, ctr, err] = Vumps(O, vumps_opts.chi_max, [], opts);
    end

end
