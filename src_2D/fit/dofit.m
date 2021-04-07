clear all; format long; close all

%% load data
names = {
    'Ising2D_g=2.5000e+00_chi=8_02_April_2021_09:38';
    'Ising2D_g=2.5000e+00_chi=11_02_April_2021_09:38';
    'Ising2D_g=2.5000e+00_chi=16_02_April_2021_09:38';
    'Ising2D_g=2.5000e+00_chi=23_02_April_2021_09:38';
    'Ising2D_g=2.5000e+00_chi=32_02_April_2021_09:38';
    'Ising2D_g=2.5000e+00_chi=45_02_April_2021_09:38';
};




delta = [];
M=[];
S=[];
T=[];
xi=[];


for i = 1:numel(names)
    data = fetch_matfiles(names{i}, struct);
    data = filter_ising_results(data, struct);
    
    marek_arr =1/3*real( (data.eps_i(:, 6)- data.eps_i(:, 1))/3+...
                         (data.eps_i(:, 5)- data.eps_i(:, 2))/2+...
                         (data.eps_i(:, 4)- data.eps_i(:, 3)));
  
    delta = [delta;marek_arr];
    M=[M; data.m ];
    S=[S;data.S];
    T=[T;data.T];
    xi=[xi; real(1./data.inv_corr_length)];
end



%% set what will be fitted: Tc, exponents: nu,beta,c, orthogonal distance to scaling function
FitTc = 1;
Fitexp = 0;
Orthdist = 0;

%% initialise
%Tc = 2 / log(1 + sqrt(2)); nu = 1; beta = 1/8; c = 1/2;

Tc = 1.2736; nu = 1; beta = 1/8; c = 1/2;


%% set dof for scaling functions
Nm = 2; Nxi = 2; NS = 2;

%% the rest
param = 1e-3 * randn(1, 9 + 8 * Nm + 8 * Nxi + 8 * NS);
param(1) = Tc; param(2) = nu; param(3) = beta; param(4) = c;

fitTc = 1;
fitexp = 0;
orthdist = 0;

% iterate until restarts don't cause improvements anymore
Er = inf; Param = param; tt = tic; tic;
while true
    %try
    [param, error] = fitt(T, delta, xi, M, S, fitTc, fitexp, orthdist, param, Nm, Nxi, NS, 200);
    if error < Er
        display(['error: ', num2str(error), ' , critical point: ', num2str(param(1), 10), ' , TC: ', num2str(Tc, 10)])
        display(['nu: ', num2str(param(2), 10), ', NU: ', num2str(nu, 10), ', beta: ', num2str(param(3), 10), ', BETA: ', num2str(beta, 10)])
        if abs(param(1) - Param(1)) / abs(param(1)) > 1e-8
            toc; tic
        end
        Param = param; Er = error;
        disp(' ')
    end
    %catch
    %display('memory error')
    %end
    if toc > 60 * 60 * 3
        break
    end
    if rand(1) < 0.2
        param = Param .* (1 + 1e-4 * randn(size(param))) + 1e-5 * randn(size(param));
    else
        param = param .* (1 + 1e-4 * randn(size(param))) + 1e-5 * randn(size(param));
    end
    if ~fitTc
        param(1) = Tc;
    end
    if ~fitexp
        param(2) = nu; param(3) = beta; param(4) = c;
    end

    fitTc = FitTc;
    fitexp = Fitexp;
    orthdist = Orthdist;

    if toc(tt) > 60 * 60 * 24 * 7
        break
    end
end
