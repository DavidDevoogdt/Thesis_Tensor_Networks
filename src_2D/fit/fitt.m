function [param, resnorm] = fitt(all_data, Fitparams, initials, maxit,free_Var,fixed_params)

    options = optimoptions('lsqnonlin', 'Display', 'iter');
    options.FunctionTolerance = 1e-12;
    options.OptimalityTolerance = 1e-12;
    options.StepTolerance = 1e-12;
    options.MaxFunctionEvaluations = 1e8;
    options.MaxIterations = maxit;
    
    if Fitparams.fitTc
        options.Display = 'iter';
    else
        options.Display = 'off';
    end
    options.UseParallel = false;
    options.Algorithm = 'Levenberg-Marquardt';

    out_data = struct();
    for i = 1:3
        Z = Fitparams.names{i};
        out_data.(Z).(free_Var) = [];
        out_data.(Z).fl = [];
        out_data.(Z).fr = [];
        out_data.(Z).gl = [];
        out_data.(Z).gr = [];
    end

    teller = 0;
    
    
    [param, resnorm, residual, exitflag, output] = lsqnonlin(@(x)fun(x), initials, [], [], options);

    function res = fun(x)
       
        teller = teller + 1;

        [x0,data,Tc,nu,beta,c,ci] =  read_params(x,Fitparams,fixed_params);
        
        delta = Fitparams.Delta_fun( all_data.eps_i, ci);

        %% make scaling functions and derivatives of scaling functions
        [data.m.fit_f, data.m.fit_df, out_data.m.fl, out_data.m.fr, out_data.m.gl, out_data.m.gr] = ScalingFunction(1, 1, x0, beta / nu, data.m.sfunparam, Fitparams.m.N, out_data.m.fl, out_data.m.fr, out_data.m.gl, out_data.m.gr);
        [data.xi.fit_f, data.xi.fit_df, out_data.xi.fl, out_data.xi.fr, out_data.xi.gl, out_data.xi.gr] = ScalingFunction(-1, 1, x0, 1 / nu, data.xi.sfunparam, Fitparams.xi.N, out_data.xi.fl, out_data.xi.fr, out_data.xi.gl, out_data.xi.gr);
        [data.S.fit_f, data.S.fit_df, out_data.S.fl, out_data.S.fr, out_data.S.gl, out_data.S.gr] = ScalingFunction(-1, 1, x0, 1 / nu, data.S.sfunparam, Fitparams.S.N, out_data.S.fl, out_data.S.fr, out_data.S.gl, out_data.S.gr);


        %% make dimensionless quantities
        X_free = all_data.(free_Var) - Tc;

        res = [];
        mask = [];

        for ii = 1:3
            %Zz = m,S or xi
            Zz = Fitparams.names{ii};

            %define fit func
            if Fitparams.logfit(ii)
                data.(Zz).fit_df = @(x)data.(Zz).fit_df(x) ./ data.(Zz).fit_f(x);
                data.(Zz).fit_f = @(x)log(data.(Zz).fit_f(x));
            end

            %get finite size scaling factors, and function values
            if ~Fitparams.subleading
                omega = 0;
                phi = 0;
                s_c_1 = 0;
                s_d_1 = 0;
            else
                omega = data.(Zz).params(2);
                phi = data.(Zz).params(4);
                s_c_1 = data.(Zz).params(1);
                s_d_1 = data.(Zz).params(3);
            end

            data.(Zz).s_corr = 1 + s_c_1 * delta.^(omega / nu);
            data.(Zz).(free_Var) = X_free .* delta.^(-1 / nu) - s_d_1 * delta.^(phi / nu);

            switch Zz
                case 'm'
                    data.m.f = log(all_data.m .* delta.^(-beta / nu) ./ data.m.s_corr);
                case 'xi'
                    data.xi.f = log(all_data.xi .* delta.^(1 / nu) ./ data.xi.s_corr);
                case 'S'
                    data.S.f = 6 / c * all_data.S + log(delta) - log(data.S.s_corr);
            end

            if ~Fitparams.logfit(ii)
                data.(Zz).f = exp(data.(Zz).f);
            end

            %normalise data before fitting
            data.(Zz).rtrange = max(data.(Zz).(free_Var)) - min(data.(Zz).(free_Var));
            data.(Zz).range = max(data.(Zz).f) - min(data.(Zz).f);

            t2 = sort(data.(Zz).(free_Var));
            n = numel(data.(Zz).(free_Var));
            data.(Zz).trange = t2(round(3/4 * n)) - t2(round(1/4 * n));

            data.(Zz).fr = data.(Zz).f / data.(Zz).range;
            data.(Zz).free_Varr = data.(Zz).(free_Var) / data.(Zz).trange;

            if isempty(out_data.(Zz).(free_Var))
                out_data.(Zz).free_Varr = data.(Zz).free_Varr;
            end

            %get error for normalised data
            if Fitparams.orthdist
                [out_data.(Zz).free_Varr, data.(Zz).res] = DistanceToLine(data.(Zz).fit_f, data.(Zz).fit_df, x0, out_data.(Zz).free_Varr, data.(Zz).free_Varr, data.(Zz).fr);
            else
                data.(Zz).res = data.(Zz).fr - data.(Zz).fit_f(data.(Zz).free_Varr);
            end

            %get mask for fitting data outside of fitting bounds
            bounds = Fitparams.bounds.(Zz);
            data.(Zz).mask = (data.(Zz).f > bounds(1)) & (data.(Zz).f < bounds(2));

            if isempty(mask)
                mask = data.(Zz).mask;
            else
                mask = mask & data.(Zz).mask;
            end

        end
        %apply mask to error measure
        for ii = 1:3
            Zz = Fitparams.names{ii};
            if Fitparams.doFit(ii)

                %apply mask
                %data.(Zz).res(~mask) = 0;

                res = [res, data.(Zz).res];
            end

        end

        %% plot
        if mod(teller - 1, numel(x)) == 0
            clf

            for ii = 1:3
                
                Zz = Fitparams.names{ii};

                free_Varr_sorted = sort(data.(Zz).free_Varr);
                f_sorted = sort(data.(Zz).f);

                fact = 3;

                xlr = [free_Varr_sorted(10), free_Varr_sorted(end - 10)] * fact;
                xl = xlr * data.(Zz).trange;
                yl = [f_sorted(10), f_sorted(end - 10)] * fact;

                data.(Zz).free_Var_lin = linspace(xlr(1), xlr(2), 200)';

                subplot(1, 3, ii)
                hold off

                if ~Fitparams.logplot(ii)
                    plot(data.(Zz).(free_Var)(mask), data.(Zz).f(mask), '.b', 'MarkerSize', 4);
                    hold on
                    plot(data.(Zz).(free_Var)(~mask), data.(Zz).f(~mask), '.r', 'MarkerSize', 4);
                    plot(data.(Zz).free_Var_lin * data.(Zz).trange, (data.(Zz).range * data.(Zz).fit_f(data.m.free_Var_lin)), '-g');

                else
                    semilogy(data.(Zz).(free_Var)(mask), data.(Zz).f(mask), '.b', 'MarkerSize', 4);
                    hold on
                    semilogy(data.(Zz).(free_Var)(~mask), data.(Zz).f(~mask), '.r', 'MarkerSize', 4);
                    semilogy(data.(Zz).free_Var_lin * data.(Zz).trange, (data.(Zz).range * data.(Zz).fit_f(data.m.free_Var_lin)), '-g');

                end
                hold off

                %set(gca, 'ylim', yl);
                set(gca, 'xlim', xl);

                xlabel(Fitparams.(Zz).xlabel);
                ylabel(Fitparams.(Zz).ylabel);
                set(gca, 'fontsize', 16)
            end

            drawnow
        end
    end
end



