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
%numbering of output states for e.g. 2x3 operator
%
% 1 4
% 2 5
% 3 6

function [rho, f, lambda] = calculate_rho(vumpser, B, M, dw)
    %vumpser: mps + environment
    %M array of @tensors elements, with some open legs
    %vumps B fixed point @uniformMps for rotated unit cell
    %dw: left upper position corresponding to left upper operator of M

    if nargin < 4
        dw = {1, 1};
    end

    d0 = dw{1} - 1;
    w0 = dw{2} - 1;

    vumpsObj = copy(vumpser);

    n = vumpsObj.width;
    m = vumpsObj.depth;

    s = size(M);
    x = s(2); y = s(1);

    lru = ff(y + 1, m); %index mps lowes row M
    lrd = ff(-lru - 2, m); %index highest mps row from below
    lc = ff(x + 1, n); %index mps column end

    %circshift unit cell for vumps and B
    vumpsObj = vumpsObj.ShiftUnitCell(d0, w0);
    B = B.ShiftUnitCell(-d0, w0); %rotate in oposite direction

    
    A = vumpsObj.mps;
    B = B.Rotate180; %make ready for application from below

    %%calculate left and right fixed point for operator between mps layers
    
    method = 'env';

    switch method
        case 'vumps'

            %get left fixed point between upper and lower fixed point for column
            [lambda1, Fl] = TransferEigs(A(lru), B(lrd), [], 5, [], 'l_LL');
            %get right fixed point between upper and lower fixed point for row,
            %shift mps's

            A2 = A.ShiftUnitCell(0, -x);
            B2 = B.ShiftUnitCell(0, -x);
            [lambda2, Fr] = TransferEigs(A2(lru), B2(lrd), [], 5, [], 'r_RR');

            %svd(Fr{1}.var);%(unitary?)
            %svd(Fl{1}.var);%(unitary?)
            lambda = (lambda1 + lambda2) / 2;

            %% create left and right tensor
            %chain left and right environments, approximate environments
            G = vumpsObj.environment;

            L = Contract({G.left(1, 1)}, {[-3, -2, -1]});
            R = G.right(1, lc);

            for i = 2:y
                L = Contract({L, G.left(ff(i, m), lc)}, {[-1:-1:-i, 1], [-(i + 2), -(i + 1), 1]});
                R = Contract({R, G.right(ff(i, m), lc)}, {[-1:-1:-i, 1], [1, -(i + 1), -(i + 2)]});
            end

            %apply transfer matrices Fr/Fl and MPS C tensors
            L = Contract({L, Fl{1}}, {[-1:-1:-(y + 1), 1], [-(y + 2), 1]});
            R = Contract({R, Fr{1}, A(1).C(lc), B(lrd).C(lc).Conj}, {[2, -(2:y + 1), 1], [1, 3], [-1, 2], [-(y + 2), 3]});

        case 'env'

            %make operator with correct size: height of M and width of unit cell
            op = TensorNone.empty(y, 0);
            for i = 1:n
                for j = 1:y
                    op(j:m:end, i) = vumpsObj.operator.mpo(j, i);
                end
            end
            op = TransferMpo(op);

            env = op.Environment(A(1), B(lrd));
            
            %reorde legs and apply Ac
            L = Contract(  { env.left(1) }, {[-(y+2):-1]});    
            R = Contract(  {env.right(lc), A(1).C(lc), B(lrd).C(lc).Conj  }, {  [1, -2:-1:-(y+1),2], [-1, 1],[-(y+2),2]   }  ); 
            
        otherwise
            error('unknown')
    end

    %% create sandwich layers and contract with left environment
    for i = 1:x

        curr = (L.legs - (y + 2));
        ext = curr + y * 2; %make room for extra legs

        ii = ff(i, m);

        tensors = cell(1, y + 2);
        tensors{1} = A(1).AL(ii);
        tensors{y + 2} = B(lrd).AL(ii).Conj;

        ord = cell(1, y + 2);

        ord{1} = [1, y + 3, -(ext + 1)];
        ord{y + 2} = [y + 2, y + 2 + (y + 1), -(ext + y + 2)];

        for j = 1:y
            ord{j + 1} = [(y + 2) + j, -(j + 1 + ext), (y + 2) + j + 1, j + 1, - (curr + 2 * j - 1), - (curr + 2 * j)];
            tensors{j + 1} = M(j, i);
        end

        L = Contract([{L}, tensors], [{[-(1:curr), 1:(y + 2)]}, ord]);
    end

    %% contract right environment

    curr = (L.legs - (y + 2));

    rho2 = Contract({L, R}, {[-(1:curr), 1:y + 2], [1:y + 2]});
    rho2 = rho2.Permute([1:2:curr, 2:2:curr]);

    %% reshape and provide handle to create operator matrix
    dd = numel(rho2.dims) / 2;
    ord = cell(1, dd);
    for i = 1:dd
        ord{i} = [-i, - (i + dd)];
    end

    tdim = prod(rho2.dims(1:dd));

    f = @(x) reshape(ncon(repmat({x}, 1, dd), ord), tdim, tdim);

    rho = reshape(rho2.var, tdim, tdim);
    %% normalise
    rho = rho / trace(rho);
end

function a = ff(a, n)
    a = mod(a - 1, n) + 1;
end
