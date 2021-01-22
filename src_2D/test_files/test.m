% function test
%     %test_H_exp
%     test_PEPO
%     %test3
% end

%function test_PEPO

function test

    disp('dd')
    fprintf("\n")

    d = 2;

    %hamiltonian setup
    S_x = [0, 1; 1, 0];
    S_y = [0, -1i; 1i, 0];
    S_z = [1, 0; 0, -1];
    I_tensor = eye(2);

    J = 1;
    g = 0.5
    %
    H_1_tensor = -J * g * S_x;
    H_2_tensor = -J * (reshape(ncon({S_z, S_z}, {[-1, -3], [-2, -4]}), [d, d, d, d]));
    %
    %
    %     H_2_tensor = -ncon({S_x, S_x}, {[-1, -3], [-2, -4]}) ...
    %             -ncon({S_y, S_y}, {[-1, -3], [-2, -4]}) ...
    %             -ncon({S_z, S_z}, {[-1, -3], [-2, -4]});
    %     H_1_tensor = zeros(d);

    %
    %     H_2_tensor = -ncon({S_x, S_x}, {[-1, -3], [-2, -4]}) ...
    %             -ncon({S_y, S_y}, {[-1, -3], [-2, -4]}) ...
    %             -ncon({S_z, S_z}, {[-1, -3], [-2, -4]});
    %    H_1_tensor = zeros(d);

    opts.testing = 0;
    opts.visualise = 0;

    pepo_order = 4;

    %T = 10.^( -2:0.5:5 )   ;

    beta_arr = 10.^(-2:1:2);
    %beta_arr=1./T;

    beta_len = size(beta_arr, 2);
    err_arr = zeros(beta_len, 1);

    for i = 1:beta_len
        beta = beta_arr(i);

        pepo = PEPO(d, -beta * H_1_tensor, ...
            -beta * H_2_tensor, ...
            pepo_order, @make_PEPO_1D, opts);

        [err, prefact] = calculate_error(pepo, 1:10, struct("numbered", true, "h_cyclic", 0));

        fprintf(" beta %.4e err %.4e \n", beta, err);
        err_arr(i) = abs(err);

        %fprintf(" beta %.4e rel err %.4e abs err %.4e \n",beta,abs(err), abs(err)* prefact );
        %err_arr(i) = abs(err);

    end

    %%
    %hold off
    figure();
    x_width = 15;
    y_width = 10;
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, x_width, y_width], 'PaperSize', [x_width, y_width])

    loglog(beta_arr, err_arr);
    %plot(T,err_arr );

    %hold on

    title(sprintf("1D transverse Ising, g=%.4f ", g));
    xlabel("$\frac{J}{k T}$", "Interpreter", "Latex");
    ylabel("err");
    ylim([0, 10])

    figure(gcf)

    %end

end
