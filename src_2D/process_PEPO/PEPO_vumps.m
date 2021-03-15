function [A, G1, lambda1, ctr, err] = PEPO_vumps(obj, chimax, maxit, name)

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

    opts.save = 1;
    opts.name = name;

    %opts.method = 'qr';

    opts.plot = 'on';
    opts.maxit = maxit;
    opts.tolfixed = 1e-10;

    %put into vumps format
    T = obj.PEPO_matrix;
    M = ncon({T}, {[1, 1, -1, -2, -3, -4]});

    o.legs = 4;
    o.group = 'none';
    o.dims = size(M);
    o.var = M;

    O.type = 'mpo';
    O.mpo = o;

    [A, G1, lambda1, ctr, err] = Vumps(O, chimax, [], opts);
end
