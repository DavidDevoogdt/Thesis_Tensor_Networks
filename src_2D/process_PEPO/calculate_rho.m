% Index convention
%
% __         __     1 6            __        __
%|  |-3   1-|  |    |/            |Fl|-1  1-|Fr|
%|Gl|-2   2-|Gr|  4-m-2   1-Al-3  |__|-2  2-|__|
%|__|-1   3-|__|   /|       |
%                 5 3       2
%
%  1- Al - n+1
%     |(2n+1)
%     |/
%  2- m  - n+2
%    /|
%2n+2 |
%    ...
%     |
%     |
%  n- Bl - 2*n
%
%

function [rho, f] = calculate_rho(A, B, G, m)

    [Al, Ar, Ca, Ac] = A{:};
    [Bl, Br, Cb, Bc] = B{:};
    [Gl, Gr] = G{:};

    %overlap between A and B
    xtens.group = 'none';
    xtens.legs = 2;
    xtens.var = rand(chi);
    xtens.dims = size(xtens.var);

    opts = struct('krylovdim', 100, 'tol', 10e-14, 'level', 1);

    [Fl, ~] = TensorEigs(@(x) TensorContract({x, Al, Bl}, {[1, 2], [1, 3, -1], [2, 3, -2]}), xtens, 1, 'lm', opts);
    [Fr, ~] = TensorEigs(@(x) TensorContract({x, Ar, Br}, {[1, 2], [-1, 3, 1], [-2, 3, 2]}), xtens, 1, 'lm', opts);

    %
    s = size(m);

    x = s(2); y = s(1);

    %create left and right env
    switch y
        case 1
            L = TensorContract({Gl, Fl}, {[1, -2, -1], [1, -3]});
            R = TensorContract({Gr, Fr, Ca, Cb}, {[2, -2, 1], [1, 3], [-1, 2], [-3, 3]});
        case 2
            L = TensorContract({Gl, Gl, Fl}, {[1, -2, -1], [2, -3, 1], [2, -4]});
            R = TensorContract({Gr, Gr, Fr, Ca, Cb}, {[3, -2, 1], [1, -3, 2], [2, 4], [-1, 3], [-4, 4]});
        otherwise
            error('generalise this')
    end

    %create sandwich layers
    for i = 1:x

        tensors = cell(1, y + 2);
        tensors{1} = Al;
        tensors{y + 2} = Bl;

        ord = cell(1, y + 2);

        ord{1} = [-1, 1, -(y + 3)];
        ord{y + 2} = [-(y + 2), y + 1, -(2 * (y + 2))];

        for j = 1:y
            ord{j + 1} = [j, -(1 + j + (y + 2)), j + 1, -(1 + j), - (2 * (y + 2) + 2 * j - 1), - (2 * (y + 2) + 2 * j)];
            tensors{j + 1} = m(j, i);
        end

        S{i} = TensorContract(tensors, ord);
    end

    %connect L, all S's and R together

    tensors = {L, S{:}, R};
    ord = cell(1, x + 2);

    ord{1} = [1:y + 2];
    ord{x + 2} = [x * (y + 2) + 1:(x + 1) * (y + 2)];

    for i = 1:x
        ord{i + 1} = [(i - 1) * (y + 2) + 1:i * (y + 2), i * (y + 2) + 1:(i + 1) * (y + 2), - (i - 1) * y, -(i + 1) * y - 1];
    end

    [rho2, ~] = TensorContract(tensors, ord);

    dd = numel(rho2.dims) / 2;
    ord = cell(1, dd);
    for i = 1:dd
        ord{i} = [-i, - (i + dd)];
    end

    tdim = prod(rho2.dims(1:dd));

    f = @(x) reshape(ncon(repmat({x}, 1, dd), ord), tdim, tdim);

    rho = reshape(rho2.var, tdim, tdim);
    %normalise
    rho = rho / trace(rho);
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

%______

%         [depth, width] = size(A{1});
%
%         AL = A{1}; AR = A{2};
%         [depth, width] = size(AL);
%         cAL = cell(depth, width); cAR = cell(depth, width);
%         for d = 1:depth
%             for w = 1:width
%                 cAL{d, w} = TensorConj(AL{d, w});
%                 cAR{d, w} = TensorConj(AR{d, w});
%             end
%         end
%
%         G = [G0{1}; G0{2}];
%         lambda = zeros(2 * depth, 1);
%         dpar = 1;
%
%         d = dpar; GLd = G(dpar, :);
%         % dd=mod(d,depth)+1;
%
%         left = MVumpsGL_ApplyFunction1(GLd{1}, d, AL, cAL, m2);
%         right = GLd{2};
%         cc = A{3}{dpar, 1};
%
%         TensorContract({GL{1, 1}, A{3}{1, 1}, GR{1, 2}, TensorConj(A{3}{2, 1})}, ...
%             {[1, 4, 2], [2, 3], [3, 4, 5], [1, 5]})
%
%         lambda = TensorContract({left, A{3}{1, 1}, right, TensorConj(A{3}{2, 1})}, ...
%             {[1, 4, 2], [2, 3], [3, 4, 5], [1, 5]});
%
%         [GLd2, ll] = TensorEigs(@(x)MVumpsGL_ApplyFunction1(x, d, AL, cAL, m), GLd{1}, 1, 'lm', opts);
%         ll = ll^(1 / width);
%
%         lambda(dpar) = ll;
%
%         %G={G(1:depth,:) G(depth+1:2*depth,:)};
%         %lambda=[prod(lambda(1:depth))^(1/depth) prod(lambda(depth+1:2*depth))^(1/depth)];
%
%         %
%         left = TensorContract({G0{1}{1, 1}}, {[-3, -2, -1]});
%         right = TensorContract({G0{2}{1, 2}, A{3}{1, 2}, TensorConj(A{3}{2, 2})}, {[1, -2, 2], [-1, 1], [-3, 2]});
%
%         row1 = TensorContract({left, A{1}{1, 1}, m, TensorConj(A{1}{2, 1})}, {[1, 2, 3], [1, 4, -1], [4, -2, 5, 2], [3, 5, -3]})
%         row2 = TensorContract({row1, A{1}{1, 2}, m, TensorConj(A{1}{2, 2})}, {[1, 2, 3], [1, 4, -1], [4, -2, 5, 2], [3, 5, -3]})
%
%         lambda = TensorContract({row2, right}, {[1, 2, 3], [1, 2, 3]})
%         %
%         %
%         %
%         %         left = TensorContract( { G0{1}{1,1}} ,{ [-3,-2,-1]  }  );
%         %         right =  TensorContract( { G0{2}{1,2} , A{3}{1,1}, TensorConj(A{3}{1,1}) } ,{ [1,-2,2] ,[-1,1],[-3,2] }  );
%         %
%         %         row1 = TensorContract( { left, A{1}{1,1}, m2 , TensorConj( A{1}{2,1})} ,{ [1,2,3], [1,4,-1] , [4,-2,5,2],[3,5,-3] }  )
%         %         row2 = TensorContract( { row1, A{1}{1,2}, m2 , TensorConj( A{1}{2,2})} ,{ [1,2,3], [1,4,-1] , [4,-2,5,2],[3,5,-3] }  )
%         %
%         %         lambda = TensorContract( { row2, right }, {[1,2,3],[1,2,3]}  )
%         %
%         %
%         %
%         %
%         left = TensorContract({G0{1}{1, 1}, G0{1}{2, 1}}, {[1, -2, -1], [-4, -3, 1]});
%         right = TensorContract({G0{2}{1, 1}, G0{2}{2, 1}}, {[-1, -2, 1], [1, -3, -4]});
%
%         %right =  TensorContract( { G0{2}{1,1}, G0{2}{2,1} , A{3}{1,2}, TensorConj(A{3}{2,1}) } ,{ [1,-2,2] , [2,-3,3],[-1,1],[-4,3] }  );
%
%         row1 = TensorContract({left, A{1}{1, 1}, m, m, TensorConj(A{1}{1, 1})}, {[1, 2, 3, 4], [1, 5, -1], [5, -2, 6, 2], [6, -3, 7, 3], [4, 7, -4]})
%         row2 = TensorContract({row1, A{4}{1, 2}, m, m, TensorConj(A{4}{1, 2})}, {[1, 2, 3, 4], [1, 5, -1], [5, -2, 6, 2], [6, -3, 7, 3], [4, 7, -4]})
%
%         lambda = TensorContract({row2, right}, {[1, 2, 3, 4], [1, 2, 3, 4]})
%         %
%         %
%         %
%         %         left = TensorContract( { G0{1}{1,1}, G0{1}{2,1} } ,{ [1,-2,-1] , [-4,-3,1] }  );
%         %         right =  TensorContract( { G0{2}{1,2}, G0{2}{2,2} , A{3}{1,1}, TensorConj(A{3}{1,1}) } ,{ [1,-2,2] , [2,-3,3],[-1,1],[-4,3] }  );
%         %
%         %         row1 = TensorContract( { left, A{1}{1,1}, m2 , m2 ,TensorConj( A{1}{1,1})} ,{ [1,2,3,4], [1,5,-1] , [5,-2,6,2],[6,-3,7,3],[4,7,-4] }  )
%         %
%         %
%         %                 www = mod(w+1,width)+1;
%         %
%         %                 d3 = mod(d-2, depth)+1;
%         %
%         %
%         %                 GL = G0{1}{d,w};
%         %                 GR = G0{2}{d,www};
%         %                 Ac = A{4}{d,ww};
%         %
%         %                 C = A{3}{d,ww};
%         %
%         %                 [B, ~, ~] = TensorEigs(@(x) TensorContract({x, GL, m, GR},{[1, 2, 5], [1, 3, -1], [-2, 4, 2, 3], [-3, 4, 5]}), TensorConj(Ac), 1, 'lm', opts);
%         %
%         %                 [B, ~, ~] = TensorEigs(@(x) get_B(x), TensorConj(Ac), 1, 'lm', opts);
%         %                 %B= TensorConj(  A{4}{d3,w}  );
%         %
%         %                 [rho, ~] = TensorContract({B, GL, Ac, o2, GR}, ...
%         %                     {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3, -1, -2], [8, 7, 6]});
%         %                 disp(rho.var)
%         %
%         %                 mag = trace(rho.var * X) / trace(rho.var)
%         %             end
%         %         end

%_______

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
%     end   %         row2 = TensorContract( { row1, A{1}{1,2}, m2 , m2 ,TensorConj(A{1}{1,2}) } ,{ [1,2,3,4], [1,5,-1] , [5,-2,6,2],[6,-3,7,3],[4,7,-4] }  )
%         %
%         %         lambda = TensorContract( { row2, right }, {[1,2,3,4],[1,2,3,4]}  )
%         %
%         %         for d=1:2
%         %             for w=1:2
%         %
%         %                 dd = mod(d, depth)+1;
%         %                 ww = mod(w, width)+1;
%         %
%         %
%         %                 www = mod(w+1,width)+1;
%         %
%         %                 d3 = mod(d-2, depth)+1;
%         %
%         %
%         %                 GL = G0{1}{d,w};
%         %                 GR = G0{2}{d,www};
%         %                 Ac = A{4}{d,ww};
%         %
%         %                 C = A{3}{d,ww};
%         %
%         %                 [B, ~, ~] = TensorEigs(@(x) TensorContract({x, GL, m, GR},{[1, 2, 5], [1, 3, -1], [-2, 4, 2, 3], [-3, 4, 5]}), TensorConj(Ac), 1, 'lm', opts);
%         %
%         %                 [B, ~, ~] = TensorEigs(@(x) get_B(x), TensorConj(Ac), 1, 'lm', opts);
%         %                 %B= TensorConj(  A{4}{d3,w}  );
%         %
%         %                 [rho, ~] = TensorContract({B, GL, Ac, o2, GR}, ...
%         %                     {[1, 2, 6], [1, 3, 4], [4, 5, 8], [5, 7, 2, 3, -1, -2], [8, 7, 6]});
%         %                 disp(rho.var)
%         %
%         %                 mag = trace(rho.var * X) / trace(rho.var)
%         %             end
%         %         end

%_______

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
