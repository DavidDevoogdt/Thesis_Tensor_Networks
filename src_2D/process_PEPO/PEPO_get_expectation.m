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
    M = ncon({T}, {[1, 1, -1, -2, -3, -4]});
    m.legs = 4;
    m.group = 'none';
    m.dims = size(M);
    m.var = M;

    M2 = ncon({T, X}, {[1, 2, -1, -2, -3, -4], [2, 1]});
    m2.legs = 4;
    m2.group = 'none';
    m2.dims = size(M2);
    m2.var = M2;
    %same but open
    O2 = ncon({T}, {[-5, -6, -1, -2, -3, -4]});
    o2.legs = 6;
    o2.group = 'none';
    o2.dims = size(O2);
    o2.var = O2;

    % calculate eigentensor under

    %get Ac equivalent from under
    opts = [];
    opts.krylovdim = 100; opts.tol = 1e-14;
    opts.level = 1;

    function a = overlap(G, A)
        a = TensorContract({G{1}, A{3}, G{2}, TensorConj(A{3})}, ...
            {[1, 4:A{1}.legs + 1, 2], [2, 3], [3, 4:A{1}.legs + 2], [1, A{1}.legs + 2]});
    end
    %
    %     function x = mpo_transf(x,AL,cAL,Ox)
    %         x = TensorContract({cAL,x,AL,Ox},{[1,2,-1],[1,3,4],[4,5,-3],[5,-2,2,3]});
    %     end
    %
    %
    %     function x = mpo_multi_tranf(x,d,AL,cAL,Ox)
    %         dd=mod(d,depth)+1;
    %
    %         [~,width]=size(AL);
    %         for ww=1:width
    %             x=TensorContract({cAL{dd,ww},x,AL{d,ww},Ox },...
    %             {[1,2,-1],[1,3,4],[4,5,-3],[5,-2,2,3]});
    %         end
    %     end
    %
    %
    %     if vumps_opts.cell_size == 1
    %         AL=A{1}; cAL=TensorConj(AL);
    %         [~,lambda]=TensorEigs(@(x) mpo_transf(x,AL,cAL,m),G0{1},1,'lm',opts);
    %
    %         G01 = G0;
    %         A1 = A;
    %
    %     else
    %
    %
    %
    %         AL=A{1};
    %         [depth,width]=size(AL);
    %         cAL=cell(depth,width);
    %         for d=1:depth
    %             for w=1:width
    %             cAL{d,w}=TensorConj(AL{d,w});
    %             end
    %         end
    %
    %         G01 = cellfun( @(x) x(1,1) ,G0  );
    %         A1 = cellfun( @(x) x(1,1) ,A  );
    %
    %         [~,lambda]=TensorEigs(@(x)mpo_multi_tranf(x,d,AL,cAL,m),  G0{1}{1} ,1,'lm',opts);
    %
    %     end

    function x = MVumpsGL_ApplyFunction1(x, d, AL, cAL, O)
        dd = mod(d, depth) + 1;
        [~, width] = size(AL);
        for ww = 1:width
            x = TensorContract({cAL{dd, ww}, x, AL{d, ww}, O}, ...
                {[1, 2, -1], [1, 3, 4], [4, 5, -3], [5, -2, 2, 3]});
        end
    end

    if vumps_opts.cell_size == 1
        GL = G0{1};
        GR = G0{2};
        Ac = A{4};
        C = A{3};

        %B = TensorConj(Ac); % works in symmetric case

        [B, ~, ~] = TensorEigs(@(x) get_Bc(x, GL, GR), TensorConj(Ac), 1, 'lm', opts);

        [rho, ~] = TensorContract({B, GL, Ac, o2, GR}, ...
            {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3, -1, -2], [8, 7, 6]});

        %mag = trace(rho.var * X) / trace(rho.var);

    else
        [depth, width] = size(A{1});

        AL = A{1}; AR = A{2};
        [depth, width] = size(AL);
        cAL = cell(depth, width); cAR = cell(depth, width);
        for d = 1:depth
            for w = 1:width
                cAL{d, w} = TensorConj(AL{d, w});
                cAR{d, w} = TensorConj(AR{d, w});
            end
        end

        G = [G0{1}; G0{2}];
        lambda = zeros(2 * depth, 1);
        dpar = 1;

        d = dpar; GLd = G(dpar, :);
        % dd=mod(d,depth)+1;

        left = MVumpsGL_ApplyFunction1(GLd{1}, d, AL, cAL, m2);
        right = GLd{2};
        cc = A{3}{dpar, 1};

        TensorContract({GL{1, 1}, A{3}{1, 1}, GR{1, 2}, TensorConj(A{3}{2, 1})}, ...
            {[1, 4, 2], [2, 3], [3, 4, 5], [1, 5]})

        lambda = TensorContract({left, A{3}{1, 1}, right, TensorConj(A{3}{2, 1})}, ...
            {[1, 4, 2], [2, 3], [3, 4, 5], [1, 5]});

        [GLd2, ll] = TensorEigs(@(x)MVumpsGL_ApplyFunction1(x, d, AL, cAL, m), GLd{1}, 1, 'lm', opts);
        ll = ll^(1 / width);

        lambda(dpar) = ll;

        %G={G(1:depth,:) G(depth+1:2*depth,:)};
        %lambda=[prod(lambda(1:depth))^(1/depth) prod(lambda(depth+1:2*depth))^(1/depth)];

        %
        left = TensorContract({G0{1}{1, 1}}, {[-3, -2, -1]});
        right = TensorContract({G0{2}{1, 2}, A{3}{1, 2}, TensorConj(A{3}{2, 2})}, {[1, -2, 2], [-1, 1], [-3, 2]});

        row1 = TensorContract({left, A{1}{1, 1}, m, TensorConj(A{1}{2, 1})}, {[1, 2, 3], [1, 4, -1], [4, -2, 5, 2], [3, 5, -3]})
        row2 = TensorContract({row1, A{1}{1, 2}, m, TensorConj(A{1}{2, 2})}, {[1, 2, 3], [1, 4, -1], [4, -2, 5, 2], [3, 5, -3]})

        lambda = TensorContract({row2, right}, {[1, 2, 3], [1, 2, 3]})
        %
        %
        %
        %         left = TensorContract( { G0{1}{1,1}} ,{ [-3,-2,-1]  }  );
        %         right =  TensorContract( { G0{2}{1,2} , A{3}{1,1}, TensorConj(A{3}{1,1}) } ,{ [1,-2,2] ,[-1,1],[-3,2] }  );
        %
        %         row1 = TensorContract( { left, A{1}{1,1}, m2 , TensorConj( A{1}{2,1})} ,{ [1,2,3], [1,4,-1] , [4,-2,5,2],[3,5,-3] }  )
        %         row2 = TensorContract( { row1, A{1}{1,2}, m2 , TensorConj( A{1}{2,2})} ,{ [1,2,3], [1,4,-1] , [4,-2,5,2],[3,5,-3] }  )
        %
        %         lambda = TensorContract( { row2, right }, {[1,2,3],[1,2,3]}  )
        %
        %
        %
        %
        left = TensorContract({G0{1}{1, 1}, G0{1}{2, 1}}, {[1, -2, -1], [-4, -3, 1]});
        right = TensorContract({G0{2}{1, 1}, G0{2}{2, 1}}, {[-1, -2, 1], [1, -3, -4]});

        %right =  TensorContract( { G0{2}{1,1}, G0{2}{2,1} , A{3}{1,2}, TensorConj(A{3}{2,1}) } ,{ [1,-2,2] , [2,-3,3],[-1,1],[-4,3] }  );

        row1 = TensorContract({left, A{1}{1, 1}, m, m, TensorConj(A{1}{1, 1})}, {[1, 2, 3, 4], [1, 5, -1], [5, -2, 6, 2], [6, -3, 7, 3], [4, 7, -4]})
        row2 = TensorContract({row1, A{4}{1, 2}, m, m, TensorConj(A{4}{1, 2})}, {[1, 2, 3, 4], [1, 5, -1], [5, -2, 6, 2], [6, -3, 7, 3], [4, 7, -4]})

        lambda = TensorContract({row2, right}, {[1, 2, 3, 4], [1, 2, 3, 4]})
        %
        %
        %
        %         left = TensorContract( { G0{1}{1,1}, G0{1}{2,1} } ,{ [1,-2,-1] , [-4,-3,1] }  );
        %         right =  TensorContract( { G0{2}{1,2}, G0{2}{2,2} , A{3}{1,1}, TensorConj(A{3}{1,1}) } ,{ [1,-2,2] , [2,-3,3],[-1,1],[-4,3] }  );
        %
        %         row1 = TensorContract( { left, A{1}{1,1}, m2 , m2 ,TensorConj( A{1}{1,1})} ,{ [1,2,3,4], [1,5,-1] , [5,-2,6,2],[6,-3,7,3],[4,7,-4] }  )
        %         row2 = TensorContract( { row1, A{1}{1,2}, m2 , m2 ,TensorConj(A{1}{1,2}) } ,{ [1,2,3,4], [1,5,-1] , [5,-2,6,2],[6,-3,7,3],[4,7,-4] }  )
        %
        %         lambda = TensorContract( { row2, right }, {[1,2,3,4],[1,2,3,4]}  )
        %
        %         for d=1:2
        %             for w=1:2
        %
        %                 dd = mod(d, depth)+1;
        %                 ww = mod(w, width)+1;
        %
        %
        %                 www = mod(w+1,width)+1;
        %
        %                 d3 = mod(d-2, depth)+1;
        %
        %
        %                 GL = G0{1}{d,w};
        %                 GR = G0{2}{d,www};
        %                 Ac = A{4}{d,ww};
        %
        %                 C = A{3}{d,ww};
        %
        %                 [B, ~, ~] = TensorEigs(@(x) TensorContract({x, GL, m, GR},{[1, 2, 5], [1, 3, -1], [-2, 4, 2, 3], [-3, 4, 5]}), TensorConj(Ac), 1, 'lm', opts);
        %
        %                 [B, ~, ~] = TensorEigs(@(x) get_B(x), TensorConj(Ac), 1, 'lm', opts);
        %                 %B= TensorConj(  A{4}{d3,w}  );
        %
        %                 [rho, ~] = TensorContract({B, GL, Ac, o2, GR}, ...
        %                     {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3, -1, -2], [8, 7, 6]});
        %                 disp(rho.var)
        %
        %                 mag = trace(rho.var * X) / trace(rho.var)
        %             end
        %         end

    end

    function x = get_Ac(x)
        x = TensorContract({GL, x, GR, m}, ...
            {[-1, 2, 1], [1, 3, 4], [4, 5, -3], [3, 5, -2, 2]});
    end

    function x = get_C(x)
        x = TensorContract({GL, x, GR}, {[-1, 3, 1], [1, 2], [2, 3, -2]});
    end

    %
    function x = get_Bc(x, GL, GR)
        x = TensorContract({x, GL, m, GR}, ...
            {[1, 2, 5], [1, 3, -1], [-2, 4, 2, 3], [-3, 4, 5]});
    end

    function x = get_C_down(x)
        x = TensorContract({x, GL, GR}, ...
            {[1, 3], [1, 2, -1], [-2, 2, 3]});
    end

    %     if p.Results.doVumps ||~isfield(save_vars, 'B')
    %         %[Ac2, ~, ~] = TensorEigs(@(x) get_Ac(x),  Ac, 1, 'lm', opts);
    %         %C2 = TensorEigs(@(x) get_C(x), C , 1, 'lm', opts);
    %
    %         %[A2,errorA2] = VumpsSolveACC(Ac2,C2, struct('method','qr','tol',1e-14)  );
    %
    %
    %         [B, ~, ~] = TensorEigs(@(x) get_Bc(x), TensorConj(Ac), 1, 'lm', opts);
    %         %C_down_2 = TensorEigs(@(x) get_C_down(x), C , 1, 'lm', opts);
    %
    %         %[B2,errorB2] = VumpsSolveACC(Bc2,C_down_2, struct('method','qr','tol',1e-14)  );
    %
    %     else
    %         B = save_vars.B;
    %     end

    %calculate <X>

    %[rho, ~] = TensorContract({B, GL, Ac, o2, GR}, ...
    %    {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3, -1, -2], [8, 7, 6]});

    mag = abs(trace(rho.var * X) / trace(rho.var));

    %get entanglement entropy

    sv = MpsSchmidtValues(C).^2;

    S = -sum(sv .* log(sv));

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
    %partion into different contributions
    %
    %     h_size= save_vars.virtual_level_sizes_horiz;
    %     v_size = save_vars.virtual_level_sizes_vert;
    %
    %     n_h_s = numel(h_size);
    %     n_v_s = numel(v_size);
    %
    %     n_v = Ac.dims(2);
    %     n_h = GL.dims(2);
    %
    %     Z_h = zeros( n_h,n_h_s,n_h);
    %     Z_v = zeros( n_v,n_v_s,n_v);
    %
    %     start_h = [0,cumsum(h_size)];
    %     start_v =  [0,cumsum(v_size)];
    %
    %     for a = 1:n_h_s
    %       Z_h( start_h(a)+1:start_h(a+1)  ,a, start_h(a)+1:start_h(a+1) ) = eye( h_size(a) );
    %     end
    %     for a = 1:n_v_s
    %       Z_v( start_v(a)+1:start_v(a+1)  ,a, start_v(a)+1:start_v(a+1) ) = eye( h_size(a) );
    %     end
    %
    %
    %     Z_virt= ncon({B.var, GL.var, Ac.var, m.var, GR.var, Z_h,Z_h,Z_v,Z_v }, ...
    %         {[3, 8, 12], [3, 2, 1], [1, 5, 9], [6, 10, 7, 4], [12, 11, 9],[2,-1,4],[10,-3,11],[5,-2,6],[7,-4,8]  },[4,10,11,9,1,2,3,5,6,7,8,12]);

    %assign output
    save_vars.A = A;
    save_vars.G0 = G0;
    save_vars.B = B;
    save_vars.PEPO_matrix = T;
    %save_vars.level_p = level_p;

    results.m = mag;

    results.S = S;
    results.chi = vumps_opts.chi_max;

    results.ftime = now;

end
