function [Xd, D] = DistanceToLine(f, df, x0, Xi, X, Y, redo)
    % determine distance d from point [x,y] to f
    % xd is the x-value of where the perpendicular line from [x,y] crosses f
    % df is the derivative of f
    % x0 is the x-value of the discontinuity of f;
    % xi is the initial guess for xd
    tolfixed = 1e-12;
    if nargin == 6
        redo = 1;
    end

    % skip the points where df is almost zero
    invg = df(Xi);
    inddo = find(abs(invg) > 1e-8);
    inddont = find(abs(invg) <= 1e-8 | isnan(invg));
    xi = Xi(inddo);
    x = X(inddo);
    y = Y(inddo);

    invg = df(xi);
    g = -1 ./ invg;
    b = y - g .* x;

    xd = getapproxcrossing(f, g, invg, b, xi, x0, tolfixed);
    er = 1;
    tel = 0;
    xdn = xd;
    ers = ones(size(xd));
    while er > tolfixed, tel = tel + 1;
        % select the unconverged points
        usexd = xd(ers > tolfixed);

        % determine slope g and offset b of proposed line that goes through point and is orhtogonal to f
        invg = df(usexd);
        g = -1 ./ invg;
        b = y(ers > tolfixed) - g .* x(ers > tolfixed);

        % determine crossing of this line with f
        xdn(ers > tolfixed) = getapproxcrossing(f, g, invg, b, usexd, x0, tolfixed);

        % update point particular errors
        ers = abs(xdn - xd);
        er = max(ers);

        % update proposed crossing point of tangent line
        xd = xdn;
        if tel > 1e2
            xd(ers > tolfixed) = x(ers > tolfixed);
            break
        end
    end

    d = sqrt((x - xd).^2 + (f(xd) - y).^2);

    % repack with the points at negligible df
    Xd(inddo) = xd;
    D(inddo) = d;
    Xd(inddont) = X(inddont);
    D(inddont) = abs(f(X(inddont)) - Y(inddont));
    if numel(Xd) > 0, Xdd = reshape(Xd, size(X)); Xd = Xdd; end
    if numel(D) > 0, D = reshape(D, size(X)); end

    % possibly consider other side of discontinuity
    indredo = find(D >= abs(X - x0));
    if redo &&~isempty(indredo)
        Xii = x0 + 1e-2 * (-1).^((Xd(indredo) - x0) > 0);
        [Xd2, D2] = DistanceToLine(f, df, x0, Xii, X(indredo), Y(indredo), 0);
        Xd(indredo) = Xd(indredo) + (Xd2 - Xd(indredo)) .* (D2 < D(indredo));
        D(indredo) = D(indredo) + (D2 - D(indredo)) .* (D2 < D(indredo));
    end
end

function xn = getapproxcrossing(f, g, invg, b, xd, x0, tolfixed)
    bt = f(xd) - invg .* xd;
    xn = (bt - b) ./ (g - invg);
    xn(xn < x0 & xd >= x0) = x0 + tolfixed;
    xn(xn > x0 & xd <= x0) = x0 - tolfixed;
end
