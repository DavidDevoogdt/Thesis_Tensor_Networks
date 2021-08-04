function [x0, data, Tc, nu, beta, c, ci] = read_params(x, Fitparams, fixed_params)
    x = real(x);

    if Fitparams.fitTc
        Tc = x(1);
    else
        Tc = TC;
    end

    if Fitparams.fitexp
        nu = x(2); beta = x(3); c = x(4);
    else
        nu = fixed_params.NU; beta = fixed_params.BETA; c = fixed_params.C;
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

    if Fitparams.dodelta == 1

        ci = x(end - 5:end);
        ci = [ci, -sum(ci)];

        ci = ci / sum(abs(ci));

    else

        ci = [];
    end

end
