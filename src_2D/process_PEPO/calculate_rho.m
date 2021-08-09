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

function [rho, f] = calculate_rho(vumpsObj, m,B,dw)

    if nargin < 4
       dw = {1,1};
    end
    
    [d0,w0] = dw{:};

    
    m = numel( vumpsObj.mps  );
    n = numel( vumpsObj.mps(1).AL ); 
    
    s = size(m);
    x = s(2); y = s(1);
    
    function a = ff(a,n)
        a = mod( a,n) + 1;
    end

    A = vumpsObj.mps;

    [lambda1, Fl] = TransferEigs(vumpsObj.mps( ff(d0,m ) ), B( ff(d0+1,m )),[], 20,[],'l_LL');
    [lambda2, Fr] = TransferEigs(vumpsObj.mps( ff(d0,m )), B( ff(d0+1,m )),[], 20,[],'r_RR');
    


    %create left and right env
    switch y
        case 1
            Gl =  vumpsObj.environment(   ff(d0,m )   ).left( ff(w0,n)  );
            Gr = vumpsObj.environment(   ff(d0+x,m )   ).right( ff(w0+x,n)  );
            
            L = Contract({ Gl   , Fl(1) }, {[1, -2, -1], [1, -3]}); 
            R = Contract({Gr, Fr, Ca, Cb}, {[2, -2, 1], [1, 3], [-1, 2], [-3, 3]});
        case 2
            L = Contract({Gl, Gl, Fl}, {[1, -2, -1], [2, -3, 1], [2, -4]});
            R = Contract({Gr, Gr, Fr, Ca, Cb}, {[3, -2, 1], [1, -3, 2], [2, 4], [-1, 3], [-4, 4]});
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



