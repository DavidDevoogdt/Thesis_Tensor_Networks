clear all; format long; close all

%% load data

[names, T_c_arr] = ising_names(2);


T_c_noise = rand * 0.01;

T_range = 0.08;

select = 1;
Tc = T_c_arr(select);
names = names{select};

all_data = struct();

%all_data.delta = [];
all_data.m = [];
all_data.S = [];
%all_data.T = [];
all_data.xi = [];
all_data.eps_i = [];

for i = 1:numel(names)
    data = fetch_matfiles(names{i}, struct);
    data = filter_ising_results(data, struct('Tbound', [Tc - T_range, Tc + T_range]));

    free_Var = data.free_var;
  
    if i == 1 
        all_data.(free_Var) = [];
    end

    all_data.eps_i = [all_data.eps_i; data.eps_i];
    all_data.m = [all_data.m; data.m];
    all_data.S = [all_data.S; data.S];
    all_data.(free_Var) = [all_data.(free_Var); data.(free_Var)];
    all_data.xi = [all_data.xi; real(1 ./ data.inv_corr_length)];
end

%% initialise
nu = 1; beta = 1/8; c = 1/2;

%universal
s_omega = 1;
s_phi = 1;
%nonuniversal
s_c = 0;
s_d = 0;

%% set what will be fitted: Tc, exponents: nu,beta,c, orthogonal distance to scaling function
bounds = struct();
bounds.S = [-Inf, Inf];
bounds.m = [-1, Inf];
bounds.xi = [-Inf, Inf];

Fitparams = struct();

Fitparams.fitTc = 1;
Fitparams.fitexp = 0;
Fitparams.orthdist = 1;
Fitparams.subleading = 0; %https://arxiv.org/pdf/cond-mat/0505194.pdf
Fitparams.doFit = [1, 1, 1];
Fitparams.names = {'m', 'xi', 'S'};
Fitparams.logplot = [0, 0, 0];
Fitparams.logfit = [0, 1, 1];
Fitparams.truncate = 1;
Fitparams.dodelta = 0;
Fitparams.Delta_fun = @eps_to_delta;

Fitparams.m.xlabel = '(T-T_c)\delta^{-1/\nu}';
Fitparams.m.ylabel = 'm \delta^{-\beta/ \nu}';

Fitparams.xi.xlabel = '(T-T_c)\delta^{-1/\nu}';
Fitparams.xi.ylabel = 'log(\xi)+log(\delta)/\nu';

Fitparams.S.xlabel = '(T-T_c)\delta^{-1/\nu}';
Fitparams.S.ylabel = ' cS/6+log(\delta)/\nu';

Fitparams.bounds = bounds;

%% set dof for scaling functioFitparams.S.N
Fitparams.m.N = 3;
Fitparams.S.N = 3;
Fitparams.xi.N = 2;

%%
%initparams
fitparams = Fitparams;
fitparams.truncate = 0;
fitparams.orthdist = 0;
fitparams.subleading = 0;
fitparams.dodelta = 0;

extra_params = 6 + 6;

%% the rest

param = 1e-2 * randn(1, 9 +12 + 8 * Fitparams.m.N + 8 * Fitparams.xi.N + 8 * Fitparams.S.N + 7);

param(1) = Tc + T_c_noise; param(2) = nu; param(3) = beta; param(4) = c;

%param(5:6)= [s_omega,s_phi];
param(5:16) = zeros(12, 1);
param(6:2:16) = 1;

%
fixed_params = struct;

if ~Fitparams.fitTc
    fixed_params.TC = param(1);
end
if ~Fitparams.fitexp
    fixed_params.NU = param(2);
    fixed_params.BETA = param(3);
    fixed_params.C = param(4);
end



ctr = 1;
% iterate until restarts don't cause improvements anymore
Er = inf; Param = param; tt = tic; tic;
while true
    %try

    if ctr == 1
        fprintf('init tries:');
        t_err = Inf;
        t_param = [];

        for i = 1:2
            param2 = Param;
            param2(17:end) = Param(17:end) .* (1 + 1e-3 * randn(size(param(17:end)))) + 1e-5 * randn(size(param(17:end)));
            %try
            [param2, error] = fitt(all_data, fitparams, param2, 50,free_Var,fixed_params);
            %catch
            %    error = t_err + 1;
            %end

            if abs(error) < t_err
                t_err = error;
                t_param = param2;
            end
        end
        param = t_param;
        error = t_err;
    else
        try
            [param, error] = fitt(all_data, fitparams, param, 20,free_Var,fixed_params);
        catch
            error = Er + 1;
        end
    end

    if error < Er
        
        [x0,data,Tcf,nuf,betaf,cf,ci] =  read_params(param,Fitparams,fixed_params);
        
        fprintf("err  %.4e \n", error);
        fprintf("Tc   %.7f (%.7f)\n", Tcf, Tc);
        fprintf("beta %.5f (%.2f)\n", betaf, beta);
        fprintf("nu   %.5f (%.1f)\n", nuf, nu);
        fprintf("c    %.5f (%.1f)\n", cf, c);
        
        if fitparams.subleading == 1
            fprintf("subleading  m:  omega %.3f  phi %.3f\n", data.m.params(2), data.m.params(4));
            fprintf("subleading  S:  omega %.3f  phi %.3f\n", data.S.params(2), data.S.params(4));
            fprintf("subleading xi:  omega %.3f  phi %.3f\n", data.xi.params(2), data.xi.params(4));
        end

        Param = param; Er = error;
        disp(' ')
    end

    param = Param .* (1 + 1e-3 * randn(size(param))) + 1e-4 * randn(size(param));

    fitparams = Fitparams;

    if toc(tt) > 60 * 60 * 24 * 7
        break
    end

    ctr = ctr + 1;
end

function delta = eps_to_delta(eps_i, c_i)

    ff= 0;

    if nargin < 2 
        ff=1;
    else
        if isempty(c_i)
            ff=1;
        end
    end

    if ff==1
       c_i = [-1, 1] / 2; 
    end
    
    delta = zeros(size(eps_i, 1), 1);

    for j = 1:numel(c_i)
        delta = delta + c_i(j) * eps_i(:, j);
    end

    delta = real(delta);

    if sum(delta) < 0
        delta = -delta;
    end
end
