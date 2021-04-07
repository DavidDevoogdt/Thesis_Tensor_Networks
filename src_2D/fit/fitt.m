function [param, resnorm] = fitt(T, delta, xi, m, S, fitTc, fitexp, orthdist, initials, Nm, Nxi, NS, maxit)
    T = T(:);
    xi = xi(:);
    m = m(:);
    delta = delta(:);
    S = S(:);

    if ~fitTc
        TC = initials(1);
    end
    if ~fitexp
        NU = initials(2);
        BETA = initials(3);
        C = initials(4);
    end

    options = optimoptions('lsqnonlin', 'Display', 'iter');
    options.FunctionTolerance = 1e-12;
    options.OptimalityTolerance = 1e-12;
    options.StepTolerance = 1e-12;
    options.MaxFunctionEvaluations = 1e8;
    options.MaxIterations = maxit;
    if fitTc
        options.Display = 'off';
    else
        options.Display = 'off';
    end
    options.UseParallel = false;
    options.Algorithm = 'Levenberg-Marquardt';
    %options.Algorithm='trust-region-reflective';

    ttxi = []; ttm = []; ttS = []; fl1 = []; fr1 = []; gl1 = []; gr1 = []; fl2 = []; fr2 = []; gl2 = []; gr2 = []; fl3 = []; fr3 = []; gl3 = []; gr3 = [];

    teller = 0; fun(initials);
    [param, resnorm, residual, exitflag, output] = lsqnonlin(@(x)fun(x), initials, [], [], options);
    param = real(param);

    if ~fitTc
        param(1) = TC;
    end
    if ~fitexp
        param(2) = NU;
        param(3) = BETA;
        param(4) = C;
    end

    teller = numel(param);
    fun(param);

    function res = fun(x)
        x = real(x);
        teller = teller + 1;

        if fitTc
            Tc = x(1);
        else
            Tc = TC;
        end

        if fitexp
            nu = x(2); beta = x(3); c = x(4);
        else
            nu = NU; beta = BETA; c = C;
        end

        %% x0 is the position of of the discontinuity of scaling functions
        x0 = x(5);

        %% make scaling functions and derivatives of scaling functions
        sfunparam = x(6:9 + 8 * Nm);
        [mf, dmf, fl1, fr1, gl1, gr1] = ScalingFunction(1, 1, x0, beta, sfunparam(1), sfunparam(2), sfunparam(3), sfunparam(4), Nm, sfunparam(5:(4 + 8 * Nm)), fl1, fr1, gl1, gr1);

        sfunparam = x(6 + 8 * Nm:9 + 8 * Nm + 8 * Nxi);
        [xif, dxif, fl2, fr2, gl2, gr2] = ScalingFunction(-1, 1, x0, nu, sfunparam(1), sfunparam(2), sfunparam(3), sfunparam(4), Nm, sfunparam(5:(4 + 8 * Nm)), fl2, fr2, gl2, gr2);

        sfunparam = x(6 + 8 * Nm + 8 * Nxi:9 + 8 * Nm + 8 * Nxi + 8 * NS);
        [Sf, dSf, fl3, fr3, gl3, gr3] = ScalingFunction(-1, 1, x0, nu, sfunparam(1), sfunparam(2), sfunparam(3), sfunparam(4), Nm, sfunparam(5:(4 + 8 * Nm)), fl3, fr3, gl3, gr3);
        dSf = @(x)dSf(x) ./ Sf(x);
        Sf = @(x)log(Sf(x));

        %% make dimensionless quantities
        t = T - Tc;

        tt = t .* delta.^(-1 / nu); trange = max(abs(tt));
        mt = (m .* delta.^(-beta / nu)); mrange = max(mt) - min(mt);
        xit = (xi .* delta.^(1 / nu)); xirange = max(xit) - min(xit);
        St = S + c / 6 * log(delta); Srange = max(St) - min(St);

        tt = tt / trange; mt = mt / mrange; xit = xit / xirange; St = St / Srange;

        if isempty(ttxi)
            ttxi = tt; ttm = tt; ttS = tt;
        end

        if orthdist
            [ttm, resm] = DistanceToLine(mf, dmf, x0, ttm, tt, mt);
            [ttxi, resxi] = DistanceToLine(xif, dxif, x0, ttxi, tt, xit);
            [ttS, resS] = DistanceToLine(Sf, dSf, x0, ttS, tt, St);
        else
            resm = mt - mf(tt);
            resxi = xit - xif(tt);
            resS = St - Sf(tt);
        end

        %% distance from points to scaling function vector
        res = [resm; resxi; resS];

        %% plot
        if mod(teller - 1, numel(x)) == 0
            ttt = linspace(min(tt), max(tt), 200)';
            tt = tt * trange; mt = mt * mrange; xit = real(xit * xirange); St = St * Srange;
            clf
            subplot(1, 3, 1)
            hold off
            plot(tt, mt, '.b', 'MarkerSize', 4);
            hold on
            plot(ttt * trange, (mrange * mf(ttt)), '-g');
            hold off
            %title('correlation length scaling function');
            xlabel('(T-T_c)\delta^{-1/\nu}');
            ylabel('m \delta^{-\beta/ \nu}');
            set(gca, 'fontsize', 16)

            subplot(1, 3, 2)
            hold off
            plot(tt, xit, '.b', 'MarkerSize', 4);
            hold on
            plot(ttt * trange, real((xirange * xif(ttt))), '-g');
            hold off
            %title('correlation length scaling function');
            xlabel('(T-T_c)\delta^{-1/\nu}');
            ylabel('\xi \delta^{1/ \nu}');
            set(gca, 'fontsize', 16)

            subplot(1, 3, 3)
            hold off
            plot(tt, St, '.b', 'MarkerSize', 4);
            hold on
            plot(ttt * trange, real((Srange * Sf(ttt))), '-g');
            hold off
            %title('correlation length scaling function');
            xlabel('(T-T_c)\delta^{-1/\nu}');
            ylabel('S + c/6 log(\delta)');
            set(gca, 'fontsize', 16)

            drawnow
        end
    end
end
