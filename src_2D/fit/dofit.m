clear all; format long; close all

%% load data

names = {{
    'Ising2D_g=0.0000e+00_chi=8_08_April_2021_09:52';
    'Ising2D_g=0.0000e+00_chi=11_08_April_2021_10:11';
    'Ising2D_g=0.0000e+00_chi=16_08_April_2021_10:43';
    'Ising2D_g=0.0000e+00_chi=23_08_April_2021_11:50';
    'Ising2D_g=0.0000e+00_chi=32_08_April_2021_14:41';
    'Ising2D_g=0.0000e+00_chi=45_08_April_2021_18:59';
    'Ising2D_g=0.0000e+00_chi=64_09_April_2021_05:59';
    'Ising2D_g=0.0000e+00_chi=91_11_April_2021_18:25';
    'Ising2D_g=0.0000e+00_chi=128_16_April_2021_16:04';
    }, {
    'Ising2D_g=1.5000e+00_chi=8_08_April_2021_10:15';
    'Ising2D_g=1.5000e+00_chi=11_08_April_2021_10:28';
    'Ising2D_g=1.5000e+00_chi=16_08_April_2021_10:48';
    'Ising2D_g=1.5000e+00_chi=23_08_April_2021_11:29';
    'Ising2D_g=1.5000e+00_chi=32_08_April_2021_12:19';
    'Ising2D_g=1.5000e+00_chi=45_08_April_2021_16:50';
    'Ising2D_g=1.5000e+00_chi=64_09_April_2021_04:35';
    'Ising2D_g=1.5000e+00_chi=91_10_April_2021_20:23';
    'Ising2D_g=1.5000e+00_chi=128_13_April_2021_18:47';
    }, {
    'Ising2D_g=2.5000e+00_chi=8_08_April_2021_09:52';
    'Ising2D_g=2.5000e+00_chi=11_08_April_2021_10:02';
    'Ising2D_g=2.5000e+00_chi=16_08_April_2021_10:17';
    'Ising2D_g=2.5000e+00_chi=23_08_April_2021_10:52';
    'Ising2D_g=2.5000e+00_chi=32_08_April_2021_12:39';
    'Ising2D_g=2.5000e+00_chi=45_08_April_2021_17:46';
    'Ising2D_g=2.5000e+00_chi=64_09_April_2021_05:52';
    'Ising2D_g=2.5000e+00_chi=91_10_April_2021_23:23';
    'Ising2D_g=2.5000e+00_chi=128_14_April_2021_21:21';
    }, {
    'Ising2D_g=2.9000e+00_chi=8_08_April_2021_09:52';
    'Ising2D_g=2.9000e+00_chi=11_08_April_2021_10:03';
    'Ising2D_g=2.9000e+00_chi=16_08_April_2021_10:06';
    }};

T_c_arr = [2 / log(1 + sqrt(2)), 1.980, 1.2737, 0.9];

T_c_noise = rand * 0.01;

T_range = 0.08;


select = 3;
Tc = T_c_arr(select);
names = names{select};

all_data = struct();

%all_data.delta = [];
all_data.m = [];
all_data.S = [];
all_data.T = [];
all_data.xi = [];
all_data.eps_i = [];



for i = 1:numel(names)
    data = fetch_matfiles(names{i}, struct);
    data = filter_ising_results(data, struct( 'Tbound', [Tc-T_range, Tc+T_range ] ));

  
    all_data.eps_i = [all_data.eps_i; data.eps_i];
    all_data.m = [all_data.m; data.m];
    all_data.S = [all_data.S; data.S];
    all_data.T = [all_data.T; data.T];
    all_data.xi = [all_data.xi; real(1 ./ data.inv_corr_length)];
end

%% initialise
%Tc = 2 / log(1 + sqrt(2)); nu = 1; beta = 1/8; c = 1/2;

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
Fitparams.subleading = 1; %https://arxiv.org/pdf/cond-mat/0505194.pdf
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
            [param2, error] = fitt(all_data, fitparams, param2, 50);
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
            [param, error] = fitt(all_data, fitparams, param, 20);
        catch
            error = Er + 1;
        end
    end

    if error < Er
        display(['error: ', num2str(error), ' , critical point: ', num2str(param(1), 10), ' , TC: ', num2str(Tc, 100)])
        display(['nu: ', num2str(param(2), 10), ', NU: ', num2str(nu, 10), ', beta: ', num2str(param(3), 10), ', BETA: ', num2str(beta, 10)])
        if fitparams.subleading == 1
            s_omega = param(6); s_phi = param(8);

            fprintf("subleading:  omega %.2f  phi %.2f\n", s_omega, s_phi);
        end

        Param = param; Er = error;
        disp(' ')
    end

    if toc > 60 * 60 * 3
        break
    end
    %if rand(1) < 0.2
    param = Param .* (1 + 1e-3 * randn(size(param))) + 1e-4 * randn(size(param));
    %else
    %    param = param .* (1 + 1e-4 * randn(size(param))) + 1e-5 * randn(size(param));
    %end

    %param = param .* (1 + 1e-4 * randn(size(param))) + 1e-4 * randn(size(param));

    if ~Fitparams.fitTc
        param(1) = Tc;
    end
    if ~Fitparams.fitexp
        param(2) = nu; param(3) = beta; param(4) = c;
    end

    fitparams = Fitparams;

    if toc(tt) > 60 * 60 * 24 * 7
        break
    end

    ctr = ctr + 1;
end


function delta = eps_to_delta(eps_i,c_i)

    if nargin <2
       %c_i = [-1/3,-1/2,-1,1,1/2,1/3]/(2*sum(1+1/2+1/3)); 
       c_i = [-1,1]/2; 
    end
    
    delta = zeros( size(eps_i,1),1 );
    
    for j=1:numel(c_i)
       delta = delta + c_i(j)* eps_i(: ,j) ;
    end

    delta = real(delta);
    
    if sum(delta)<0
       delta = -delta; 
    end
end

