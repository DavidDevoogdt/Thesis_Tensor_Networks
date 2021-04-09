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
Tc = 1.2736;
% 
% names = {
%     'Ising2D_g=1.5000e+00_chi=8_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=11_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=16_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=23_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=32_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=45_02_April_2021_11:17';
%     'Ising2D_g=1.5000e+00_chi=64_02_April_2021_11:17';
%     };
% Tc = 1.98;

% names ={
%     'Ising2D_g=0.0000e+00_chi=8_02_April_2021_11:46';
%     'Ising2D_g=0.0000e+00_chi=11_02_April_2021_11:46';
%     'Ising2D_g=0.0000e+00_chi=16_02_April_2021_11:46';
%     'Ising2D_g=0.0000e+00_chi=23_02_April_2021_11:46';
%     'Ising2D_g=0.0000e+00_chi=32_02_April_2021_11:46';
%     'Ising2D_g=0.0000e+00_chi=45_02_April_2021_11:46';
%     }
% Tc = 2 / log(1 + sqrt(2));


names = {{
    'Ising2D_g=0.0000e+00_chi=8_08_April_2021_09:52';
    'Ising2D_g=0.0000e+00_chi=11_08_April_2021_10:11';
    'Ising2D_g=0.0000e+00_chi=16_08_April_2021_10:43';
    'Ising2D_g=0.0000e+00_chi=23_08_April_2021_11:50';
    'Ising2D_g=0.0000e+00_chi=32_08_April_2021_14:41';
    'Ising2D_g=0.0000e+00_chi=45_08_April_2021_18:59';
    'Ising2D_g=0.0000e+00_chi=64_09_April_2021_05:59';
    },{
    'Ising2D_g=1.5000e+00_chi=8_08_April_2021_10:15';
    'Ising2D_g=1.5000e+00_chi=11_08_April_2021_10:28';
    'Ising2D_g=1.5000e+00_chi=16_08_April_2021_10:48';
    'Ising2D_g=1.5000e+00_chi=23_08_April_2021_11:29';
    'Ising2D_g=1.5000e+00_chi=32_08_April_2021_12:19';
    'Ising2D_g=1.5000e+00_chi=45_08_April_2021_16:50';
    'Ising2D_g=1.5000e+00_chi=64_09_April_2021_04:35';
    },{
    %'Ising2D_g=2.5000e+00_chi=8_08_April_2021_09:52';
    %'Ising2D_g=2.5000e+00_chi=11_08_April_2021_10:02';
    'Ising2D_g=2.5000e+00_chi=16_08_April_2021_10:17';
    'Ising2D_g=2.5000e+00_chi=23_08_April_2021_10:52';
    'Ising2D_g=2.5000e+00_chi=32_08_April_2021_12:39';
    'Ising2D_g=2.5000e+00_chi=45_08_April_2021_17:46';
    'Ising2D_g=2.5000e+00_chi=64_09_April_2021_05:52';
    },{
    'Ising2D_g=2.9000e+00_chi=8_08_April_2021_09:52';
    'Ising2D_g=2.9000e+00_chi=11_08_April_2021_10:03';
    'Ising2D_g=2.9000e+00_chi=16_08_April_2021_10:06';
}};

T_c_arr = [  2/log( 1+ sqrt(2) )   ,1.980,1.2737, 0.9  ];

T_c_noise = rand*0.01;


select = 1;
Tc = T_c_arr(select);
names = names{select};

delta = [];
M = [];
S = [];
T = [];
xi = [];

for i = 1:numel(names)
    data = fetch_matfiles(names{i}, struct);
    data = filter_ising_results(data, struct);

    marek_arr = real((data.eps_i(:, 6) - data.eps_i(:, 1)) / 3 + ...
        (data.eps_i(:, 5) - data.eps_i(:, 2)) / 2 + ...
        (data.eps_i(:, 4) - data.eps_i(:, 3)));

    delta = [delta; marek_arr];
    M = [M; data.m];
    S = [S; data.S];
    T = [T; data.T];
    xi = [xi; real(1 ./ data.inv_corr_length)];
end

%% initialise
%Tc = 2 / log(1 + sqrt(2)); nu = 1; beta = 1/8; c = 1/2;

nu = 1; beta = 1/8; c = 1/2;

%universal
s_omega = 1;
s_phi = 1;
%nonuniversal
s_c = 1e-4;
s_d = 1e-4;

%% set dof for scaling functions
Nm = 3; Nxi = 3; NS = 2;


%% set what will be fitted: Tc, exponents: nu,beta,c, orthogonal distance to scaling function
FitTc = 1;
Fitexp = 0;
Orthdist = 1;
Subleading = 1; %https://arxiv.org/pdf/cond-mat/0505194.pdf

FitS=1;

%first round
fitTc = 1;
fitexp = 0;
orthdist = 0;
subleading = 0;


extra_params = 6+6;

%% the rest
if FitS==1
    param = 1e-3 * randn(1, 9 +12+ 8 * Nm + 8 * Nxi + 8 * NS);
else
    param = 1e-3 * randn(1, 9 +12+ 8 * Nm + 8 * Nxi);
end


param(1) = Tc+T_c_noise; param(2) = nu; param(3) = beta; param(4) = c;



%param(5:6)= [s_omega,s_phi];
param(5:16) = zeros(12,1);
param(6:2:16)=1;



% iterate until restarts don't cause improvements anymore
Er = inf; Param = param; tt = tic; tic;
while true
    try
    [param, error] = fitt(T, delta, xi, M, S, fitTc, fitexp,FitS, orthdist, subleading, param, Nm, Nxi, NS, 50);
    catch
        error = Er+1;
    end
    
    if error < Er
        display(['error: ', num2str(error), ' , critical point: ', num2str(param(1), 10), ' , TC: ', num2str(Tc, 10)])
        display(['nu: ', num2str(param(2), 10), ', NU: ', num2str(nu, 10), ', beta: ', num2str(param(3), 10), ', BETA: ', num2str(beta, 10)])
        if subleading == 1
            s_omega = param(6);s_phi = param(8);
            
            fprintf("subleading:  omega %.2f  phi %.2f\n",s_omega,s_phi);
        end
        
        Param = param; Er = error;
        disp(' ')
    end


    if toc > 60 * 60 * 3
        break
    end
    %if rand(1) < 0.2
        param = Param .* (1 + 1e-3 * randn(size(param))) + 1e-5 * randn(size(param));
    %else
    %    param = param .* (1 + 1e-4 * randn(size(param))) + 1e-5 * randn(size(param));
    %end
    
    %param = param .* (1 + 1e-4 * randn(size(param))) + 1e-4 * randn(size(param));
    
    if ~fitTc
        param(1) = Tc;
    end
    if ~fitexp
        param(2) = nu; param(3) = beta; param(4) = c;
    end

    fitTc = FitTc;
    fitexp = Fitexp;
    orthdist = Orthdist;
    subleading = Subleading;

    if toc(tt) > 60 * 60 * 24 * 7
        break
    end
end
