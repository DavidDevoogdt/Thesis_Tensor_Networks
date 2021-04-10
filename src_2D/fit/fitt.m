function [param, resnorm] = fitt(all_data, Fitparams, initials, maxit)

    if ~Fitparams.fitTc
        TC = initials(1);
    end
    if ~Fitparams.fitexp
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
        out_data.(Z).temp = [];
        out_data.(Z).fl = [];
        out_data.(Z).fr = [];
        out_data.(Z).gl = [];
        out_data.(Z).gr = [];
    end

    teller = 0; fun(initials);
    [param, resnorm, residual, exitflag, output] = lsqnonlin(@(x)fun(x), initials, [], [], options);
    param = real(param);

    if ~Fitparams.fitTc
        param(1) = TC;
    end
    if ~Fitparams.fitexp
        param(2) = NU;
        param(3) = BETA;
        param(4) = C;
    end

    teller = numel(param);
    fun(param);

    function res = fun(x)
        x = real(x);
        teller = teller + 1;

        if Fitparams.fitTc
            Tc = x(1);
        else
            Tc = TC;
        end

        if Fitparams.fitexp
            nu = x(2); beta = x(3); c = x(4);
        else
            nu = NU; beta = BETA; c = C;
        end

        data = struct();

        %% read out params 
        x0 = x(5 + 12);
        
        data.m.sfunparam = x(12 + 6:12 + 9 + 8 * Fitparams.m.N);
        data.xi.sfunparam = x(12 + 6 + 8 * Fitparams.m.N:12 + 9 + 8 * Fitparams.m.N + 8 * Fitparams.xi.N);
        data.S.sfunparam = x(6 + 12 + 8 * Fitparams.m.N + 8 * Fitparams.xi.N:9 + 12 + 8 * Fitparams.m.N + 8 * Fitparams.xi.N + 8 * Fitparams.S.N);
        
        data.m.params = x(5:8);
        data.xi.params = x(9:12);
        data.S.params = x(13:16);
        
        %% make scaling functions and derivatives of scaling functions
        [data.m.fit_f, data.m.fit_df, out_data.m.fl, out_data.m.fr, out_data.m.gl, out_data.m.gr] = ScalingFunction(1, 1, x0, beta / nu, data.m.sfunparam, Fitparams.m.N, out_data.m.fl, out_data.m.fr, out_data.m.gl, out_data.m.gr);
        [data.xi.fit_f, data.xi.fit_df, out_data.xi.fl, out_data.xi.fr, out_data.xi.gl, out_data.xi.gr] = ScalingFunction(-1, 1, x0, 1 / nu, data.xi.sfunparam, Fitparams.xi.N, out_data.xi.fl, out_data.xi.fr, out_data.xi.gl, out_data.xi.gr);
        [data.S.fit_f, data.S.fit_df, out_data.S.fl, out_data.S.fr, out_data.S.gl, out_data.S.gr] = ScalingFunction(-1, 1, x0, 1 / nu, data.S.sfunparam, Fitparams.S.N, out_data.S.fl, out_data.S.fr, out_data.S.gl, out_data.S.gr);
        
        %do log

       
        %% make dimensionless quantities
        t = all_data.T - Tc;
        
        omega = x(6);
        phi = x(8);
        
        res = [];
        
        for ii = 1:3
            Zz = Fitparams.names{ii};
            
            if Fitparams.logfit(ii)
                data.(Zz).fit_df = @(x)data.(Zz).fit_df(x) ./ data.(Zz).fit_f(x);
                data.(Zz).fit_f = @(x)log(data.(Zz).fit_f(x));
            end

            if ~Fitparams.subleading
                s_c_1 =0;
                s_d_1 =0;
            else        
                s_c_1 = data.(Zz).params(1);
                s_d_1 = data.(Zz).params(3);
            end
            
            data.(Zz).s_corr = 1 + s_c_1 * all_data.delta.^(omega / nu);
            data.(Zz).temp = t .* all_data.delta.^(-1 / nu) - s_d_1 * all_data.delta.^(phi / nu);
                       

            switch Zz
                case 'm'
                    data.m.f = log(all_data.m .* all_data.delta.^(-beta / nu) ./ data.m.s_corr );
                case 'xi'
                    data.xi.f = log( all_data.xi .* all_data.delta.^(1 / nu) ./ data.xi.s_corr );
                case 'S'
                    data.S.f = 6 / c * all_data.S + log( all_data.delta) - log(data.S.s_corr) ;
            end

            if ~Fitparams.logfit(ii) 
                data.(Zz).f = exp( data.(Zz).f );
            end
                
            
            
            data.(Zz).trange = max(abs(data.(Zz).temp));
            data.(Zz).range = max(data.(Zz).f) - min(data.(Zz).f);

            data.(Zz).fr = data.(Zz).f / data.(Zz).range;
            data.(Zz).tempr = data.(Zz).temp / data.(Zz).trange;

            if isempty(out_data.(Zz).temp)
                out_data.(Zz).tempr = data.(Zz).tempr;
            end

            if Fitparams.orthdist
                [out_data.(Zz).tempr, data.(Zz).res] = DistanceToLine(data.(Zz).fit_f, data.(Zz).fit_df, x0, out_data.(Zz).tempr, data.(Zz).tempr, data.(Zz).fr);
            else
                data.(Zz).res = data.(Zz).fr - data.(Zz).fit_f(data.(Zz).tempr);
            end
      

            bounds = Fitparams.bounds.(Zz);

            data.(Zz).mask = (data.(Zz).f > bounds(1)) & (data.(Zz).f < bounds(2));

            %apply mask
            data.(Zz).res(~data.(Zz).mask) = 0;
            
        
            
            if Fitparams.doFit(ii)
                res = [res, data.(Zz).res];
            end
            
        end

           %% plot
        if mod(teller - 1, numel(x)) == 0   
            clf

            for ii = 1:3
                Zz = Fitparams.names{ii};
                data.(Zz).temp_lin = linspace(min(data.(Zz).tempr), max(data.(Zz).tempr), 200)';
                
                subplot(1, 3, ii)
                hold off
                
                if ~Fitparams.logplot(ii) 
     
                    plot(data.(Zz).temp(  data.(Zz).mask ), data.(Zz).f( data.(Zz).mask), '.b', 'MarkerSize', 4);
                    hold on
                    plot(data.(Zz).temp(  ~data.(Zz).mask ), data.(Zz).f( ~data.(Zz).mask), '.r', 'MarkerSize', 4);
                    plot(data.(Zz).temp_lin * data.(Zz).trange, (data.(Zz).range * data.(Zz).fit_f(data.m.temp_lin)), '-g');
                else
                    semilogy(data.(Zz).temp(  data.(Zz).mask ), data.(Zz).f( data.(Zz).mask), '.b', 'MarkerSize', 4);
                    hold on
                    semilogy(data.(Zz).temp(  ~data.(Zz).mask ), data.(Zz).f( ~data.(Zz).mask), '.r', 'MarkerSize', 4);
                    semilogy(data.(Zz).temp_lin * data.(Zz).trange, (data.(Zz).range * data.(Zz).fit_f(data.m.temp_lin)), '-g');
                end
                
              
                
                hold off
                
                xlabel( Fitparams.(Zz).xlabel );
                ylabel( Fitparams.(Zz).ylabel );
                set(gca, 'fontsize', 16)
            end
            
            drawnow
        end

    end
end
