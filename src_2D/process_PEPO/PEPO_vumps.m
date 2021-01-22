function [A, B, G1, lambda1] = vumps(obj, chimax)

    %todo check these params
    opts.charges = 'regular';
    opts.dynamical = 'off';
    opts.dyncharges = 0;
    opts.schmidtcut = 1e-10;
    opts.chimax = 350;
    %opts.disp='iter';
    opts.disp = 'none';
    opts.tolmax = 1e-4; %1e-4
    opts.tolfactor = 1e4;
    opts.minit = 1;
    opts.dyniter = 5;
    opts.truncate = 0;
    opts.method = 'vumps';
    opts.save = 0;

    %opts.method = 'qr';

    opts.plot = 'on';
    opts.maxit = 1000;
    opts.tolfixed = 1e-12;

    %put into vumps format

    T = obj.PEPO_matrix;

    hdim = size(T, 3);
    vdim = size(T, 4);

    %upper vumps zipper
    M = ncon({T}, {[1, 1, -1, -2, -3, -4]});
    %lower vumps zipper

    o.legs = 4;
    o.group = 'none';
    o.dims = size(M);
    o.var = M;

    O.type = 'mpo';
    O.mpo = o;

    [A, G1, lambda1, ~, ~] = Vumps(O, chimax, [], opts);

    %correct estimate for inversion sym?

    GL = G1{1}; GR = G1{2}; Ac = A{4};

    m.legs = 4;
    m.group = 'none';
    m.dims = size(M);
    m.var = M;

    opts.krylovdim = 100; opts.tol = 1e-14;
    %opts.disp='iter-detailed'; %opts.reorth='force';

    function x = vumps_under(x)
        x = TensorContract({x, GL, m, GR}, {[1, 2, 5], [1, 3, -1], [-2, 4, 2, 3], [-3, 4, 5]});
    end

    [B, lambda2, err2] = TensorEigs(@(x) vumps_under(x), TensorConj(Ac), 1, 'lm', opts);

    %[B_l,C_l,~]=TensorDecRight(Bc,'polar');

    %[B_r,C_r,~]=TensorDecLeft(Bc,'polar');

    %             accopts.method = 'qr';
    %             accopts.tol = 1e-16;
    %
    %             [B,err]=VumpsSolveACC(B2,B_c,accopts);

end
