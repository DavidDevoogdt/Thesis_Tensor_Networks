function [param, resnorm] = fitt(T, delta, xi, m, S, fitTc, fitexp,FitS, orthdist, subleading, initials, Nm, Nxi, NS, maxit)
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
    options.StepTolerance = 1e-9;
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
        x0 = x(5+12);

        %% make scaling functions and derivatives of scaling functions
        sfunparam = x(12+6:12+9 + 8 * Nm);
        [mf, dmf, fl1, fr1, gl1, gr1] = ScalingFunction(1, 1, x0, beta/nu, sfunparam(1), sfunparam(2), sfunparam(3), sfunparam(4), Nm, sfunparam(5:(4 + 8 * Nm)), fl1, fr1, gl1, gr1);

        sfunparam = x(12+6+ 8 * Nm:12+9+ 8 * Nm + 8 * Nxi);
        [xif, dxif, fl2, fr2, gl2, gr2] = ScalingFunction(-1, 1, x0, 1/nu, sfunparam(1), sfunparam(2), sfunparam(3), sfunparam(4), Nxi, sfunparam(5:(4 + 8 * Nxi)), fl2, fr2, gl2, gr2);

        if FitS == 1   
            sfunparam = x(6+12 + 8 * Nm + 8 * Nxi:9+12 + 8 * Nm + 8 * Nxi + 8 * NS);
            [Sf, dSf, fl3, fr3, gl3, gr3] = ScalingFunction(-1, 1, x0, 1/nu, sfunparam(1), sfunparam(2), sfunparam(3), sfunparam(4), NS, sfunparam(5:(4 + 8 * NS)), fl3, fr3, gl3, gr3);
            %dSf = @(x)dSf(x) ./ Sf(x);
            %Sf = @(x)log(Sf(x));
        end
        %% make dimensionless quantities
        t = T - Tc;

        if ~subleading
            tt = t .* delta.^(-1 / nu); trange = max(abs(tt));
            s_corr = 1+0*delta;
            
            tt1=tt;tt2=tt;tt3=tt;
            trange1=trange;trange2=trange;trange3=trange;
            s_corr1=s_corr;s_corr2=s_corr;s_corr3=s_corr;

        else
           
            
           
            
            pp = x(5:8); s_c_1=pp(1);omega=pp(2);s_d_1=pp(3);phi=pp(4);
            s_corr1 = 1+s_c_1*delta.^(omega/nu)  ;
            tt1 = t .* delta.^(-1 / nu) - s_d_1*delta.^(phi/nu) ; trange1 = max(abs(tt1));
            
            pp = x(9:12); s_c_1=pp(1);s_d_1=pp(3);
            s_corr2 =  1+s_c_1*delta.^(omega/nu)  ;
            tt2 = t .* delta.^(-1 / nu) - s_d_1*delta.^(phi/nu) ; trange2 = max(abs(tt2));
            
            if FitS == 1 
                pp = x(13:16);  s_c_1=pp(1);s_d_1=pp(3);
                s_corr3 = 1+s_c_1*delta.^(omega/nu)  ;
                tt3 = t .* delta.^(-1 / nu) - s_d_1*delta.^(phi/nu) ; trange3 = max(abs(tt3));
            else
                tt3 = t .* delta.^(-1 / nu)  ; trange3 = max(abs(tt3));
            end

        end

        mt = (m .* delta.^(-beta / nu))./s_corr1   ; mrange = max(mt) - min(mt);
        xit = (xi .* delta.^(1 / nu))./s_corr2; xirange = max(xit) - min(xit);
        
        if FitS == 1 
            St = exp( 6/c*S ).*delta ./s_corr3; Srange = max(St) - min(St);
            St = St / Srange;
        else
            St = exp( 6/c*S ).*delta; Srange = max(St) - min(St);
            St = St / Srange;
        end
       
        
        mt = mt / mrange; xit = xit / xirange; 

        if isempty(ttxi)
            ttxi = tt1; ttm = tt2;ttS = tt3;
           
        end

        if orthdist
            [ttm, resm] = DistanceToLine(mf, dmf, x0, ttm, tt1, mt);
            [ttxi, resxi] = DistanceToLine(xif, dxif, x0, ttxi, tt2, xit);
            if FitS == 1 
                [ttS, resS] = DistanceToLine(Sf, dSf, x0, ttS, tt3, St);
            end
        else
            resm = mt - mf(tt1);
            resxi = xit - xif(tt2);
            if FitS == 1 
                resS = St - Sf(tt3);
            end
        end

        %% distance from points to scaling function vector
        if FitS == 1 
            res = [resm; resxi; resS];
        else
            res = [resm; resxi];
        end

        %% plot
        if mod(teller - 1, numel(x)) == 0
            ttt1 = linspace(min(tt1), max(tt1), 200)';
            ttt2 = linspace(min(tt2), max(tt2), 200)';
            
            ttt3 = linspace(min(tt3), max(tt3), 200)';
            
            tt1 = tt1 * trange1;
            tt2 = tt2 * trange2;
            tt3 = tt3 * trange3;
            
            mt = mt * mrange; xit = xit * xirange; 
            
            St = St * Srange;
            
            clf
            subplot(1, 3, 1)
            hold off
            plot(tt1, mt, '.b', 'MarkerSize', 4);
            hold on
            plot(ttt1 * trange1, (mrange * mf(ttt1)), '-g');
            hold off
            %title('correlation length scaling function');
            xlabel('(T-T_c)\delta^{-1/\nu}');
            ylabel('m \delta^{-\beta/ \nu}');
            set(gca, 'fontsize', 16)

            subplot(1, 3, 2)
            hold off
            plot(tt2, xit, '.b', 'MarkerSize', 4);
            hold on
            plot(ttt2 * trange2, real((xirange * xif(ttt2))), '-g');
            hold off
            %title('correlation length scaling function');
            xlabel('(T-T_c)\delta^{-1/\nu}');
            ylabel('\xi \delta^{1/ \nu}');
            set(gca, 'fontsize', 16)

            subplot(1, 3, 3)
            hold off
            semilogy(tt3, St, '.b', 'MarkerSize', 4);
            if FitS == 1 

                hold on
                semilogy(ttt3 * trange3, real((Srange * Sf(ttt3))), '-g');
                hold off
            end
            
            %title('correlation length scaling function');
            xlabel('(T-T_c)\delta^{-1/\nu}');
            ylabel( ' e^{ c S/6}  \delta^{1/\nu}');
            set(gca, 'fontsize', 16)

            drawnow
        end
    end
end
