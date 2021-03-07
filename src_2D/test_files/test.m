% function test
%     %test_H_exp
%     test_PEPO
%     %test3
% end

%function test_PEPO

% function test
%
%     disp('dd')
%     fprintf("\n")
%
%     d = 2;
%
%     %hamiltonian setup
%     S_x = [0, 1; 1, 0];
%     S_y = [0, -1i; 1i, 0];
%     S_z = [1, 0; 0, -1];
%     I_tensor = eye(2);
%
%     %handle = @make_PEPO_1D;
%     handle = @make_PEPO_1D_double;
%
%     J = 1;
%     g = 0.5;
%     %
%     H_1_tensor = -J * g * S_x;
%     H_2_tensor = -J * (reshape(ncon({S_z, S_z}, {[-1, -3], [-2, -4]}), [d, d, d, d]));
%     %
%     %
%     %     H_2_tensor = -ncon({S_x, S_x}, {[-1, -3], [-2, -4]}) ...
%     %             -ncon({S_y, S_y}, {[-1, -3], [-2, -4]}) ...
%     %             -ncon({S_z, S_z}, {[-1, -3], [-2, -4]});
%     %     H_1_tensor = zeros(d);
%
%     %
%     %     H_2_tensor = -ncon({S_x, S_x}, {[-1, -3], [-2, -4]}) ...
%     %         -ncon({S_y, S_y}, {[-1, -3], [-2, -4]}) ...
%     %         -ncon({S_z, S_z}, {[-1, -3], [-2, -4]});
%     %     H_1_tensor = zeros(d);
%
%     opts.testing = 0;
%     opts.visualise = 0;
%
%     pepo_order = 5;
%
%     %T = 10.^( -2:0.5:5 )   ;
%
%     beta_arr = 10.^(-3:0.2:2);
%     %beta_arr=1./T;
%
%     beta_len = size(beta_arr, 2);
%     err_arr = zeros(beta_len, 1);
%
%     for i = 1:beta_len
%         beta = beta_arr(i);
%
%         pepo = PEPO(d, -beta * H_1_tensor, ...
%             -beta * H_2_tensor, ...
%             pepo_order, handle, opts);
%
%         [err, ~] = calculate_error(pepo, 1:10, struct("numbered", true, "h_cyclic", 1),1);
%
%         fprintf(" beta %.4e err %.4e \n", beta, err);
%         err_arr(i) = abs(err);
%
%         %fprintf(" beta %.4e rel err %.4e abs err %.4e \n",beta,abs(err), abs(err)* prefact );
%         %err_arr(i) = abs(err);
%
%     end
%
%     %%
%     %hold off
%     figure();
%     x_width = 15;
%     y_width = 10;
%     set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, x_width, y_width], 'PaperSize', [x_width, y_width])
%
%     loglog(beta_arr, err_arr);
%     %plot(T,err_arr );
%
%     %hold on
%
%     title(sprintf("1D transverse Ising, g=%.4f ", g));
%     xlabel("$\frac{J}{k T}$", "Interpreter", "Latex");
%     ylabel("err");
%     ylim([0, 10])
%
%     figure(gcf)
%
%     %end
%
% end

function test

    %simulatiemodellen = [ "t_ising","Heisenberg_2D",  "Heisenberg_2D_X", "random", "random"];
    simulatiemodellen = ["t_ising"];
    %simulatiemodellen = ["Heisenberg_2D",  "Heisenberg_2D_X"];
    %simulatiemodellen = [ "Heisenberg_2D","t_ising",  "Heisenberg_2D_X", "random", "random"];

    models_len = size(simulatiemodellen, 2);

    for round = 1:models_len
        %change this
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        generate_opts.testing = 0;
        generate_opts.MPO_type = "matrix";

        model_opts.g = 0.5;
        %[simul,H_1_tensor,H_2_tensor,opt4,d] = models('Heisenberg_2D',model_opts);
        %[simul,H_1_tensor,H_2_tensor,opt4,d] = models('t_ising',model_opts);
        [simul, H_1_tensor, H_2_tensor, opt4, d] = models(simulatiemodellen(round), model_opts);

        %     H_2_tensor = H_2_tensor+...
        %                    0.5* reshape( ncon( {H_1_tensor,eye(d)}, {[-1,-3],[-2,-4]}), [d,d,d,d])+...
        %                    0.5* reshape( ncon( {eye(d),H_1_tensor}, {[-1,-3],[-2,-4]}), [d,d,d,d]);
        %      H_1_tensor = 0*eye(d);

        %opt4.single_threshold = 1e-12;
        %
        %
        % %simulatiemodellen = ["random","random","random","Heisenberg_2D","t_ising","Heisenberg_2D_X","random","random"];
        % simulatiemodellen = [1e-10,1e-11,1e-12,1e-13,1e-14];
        % models_len = size(simulatiemodellen,2);
        %
        % generate_opts.testing=0;
        % generate_opts.MPO_type = "matrix";
        %
        % opt3 = struct([]);
        %
        %
        % model_opts.g =  0.5;
        % %[simul,H_1_tensor,H_2_tensor,opt4,d] = models('Heisenberg_2D',model_opts);
        % %[simul,H_1_tensor,H_2_tensor,opt4,d] = models('t_ising',model_opts);
        % [simul,H_1_tensor,H_2_tensor,opt4,d] = models( "t_ising"  ,model_opts);
        %
        %     for round = 1:models_len
        %
        %     opt4.single_threshold =simulatiemodellen(round);
        %

        simul.Order_arr = [3, 4, 5, 6];
        simul.types = [2, 3];
        simul.M = 9;
        simul.beta_arr = 10.^(-2:0.05:1.5);
        %simul.beta_arr = 10.^(  -0:0.05:1);
        simul.cyclic = 1;

        %opt5.method = "svd";

        opt2 = {};
        opt3.SparseArr = 1;
        opt5.SparseArr = 1;
        %opt5= struct([]);

        compare_opts.ref = 2;
        compare_opts.cyclic = simul.cyclic;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%plot stuff
        order_size = size(simul.Order_arr, 2);
        line_spec = ["-", "--", "-.", ":", "-", "--", "-.", ":"];
        alphas = [1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5];
        colors = {[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 1]}; %reserved for specific type

        legend_Arr = cell(size(simul.types, 2) * order_size, 1);
        legend_Arr(:) = {"todo"};
        beta_len = size(simul.beta_arr, 2);

        plot_counter = 1;
        %hold off
        figure();
        x_width = 15;
        y_width = 10;
        set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, x_width, y_width], 'PaperSize', [x_width, y_width])
        %%%

        %

        %loop over different orders
        for j = 1:order_size
            Order = simul.Order_arr(j);
            plot_structure = cell(5, beta_len);

            %loop over temps
            for i = 1:beta_len
                beta = simul.beta_arr(i);

                %mpo_base_matrix = mpo_base.H_exp(simul.M-1, 1, simul.cyclic); %buffered, not calculted again every time

                fprintf("M %d beta %.4e order %d", simul.M, beta, Order);

                %loop of simulation types
                for t = 1:size(simul.types, 2)
                    switch simul.types(t)
                        case 2
                            handle = @make_PEPO_1D;
                            pepo = PEPO(d, -beta * H_1_tensor, ...
                                -beta * H_2_tensor, ...
                                Order, handle, struct());

                            fprintf(".");
                            [err_02, ~] = calculate_error(pepo, 1:10, struct("numbered", true, "h_cyclic", simul.cyclic), 1);

                            fprintf(" err 02 %.4e", err_02);
                            plot_structure{2, i} = err_02;
                        case 3
                            handle = @make_PEPO_1D_double;
                            pepo = PEPO(d, -beta * H_1_tensor, ...
                                -beta * H_2_tensor, ...
                                Order, handle, struct());

                            fprintf(".");
                            [err_03, ~] = calculate_error(pepo, 1:10, struct("numbered", true, "h_cyclic", simul.cyclic), 1);

                            fprintf(" err 03 %.4e", err_03);
                            plot_structure{3, i} = err_03;
                        case 4
                            error("")
                        case 5
                            error("")
                        otherwise
                            error("unknown type")
                    end
                end

                fprintf("\n");

            end

            %plotting loop for current order
            for t = 1:size(simul.types, 2)
                switch simul.types(t)
                    case 2
                        colour = colors{1};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(cell2mat(plot_structure(2, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("A:%d", Order);
                        plot_counter = plot_counter + 1;
                        hold on
                    case 4
                        colour = colors{3};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(cell2mat(plot_structure(4, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("A:%d", Order);
                        plot_counter = plot_counter + 1;
                        hold on
                    case 3
                        colour = colors{2};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(cell2mat(plot_structure(3, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("C:%d", Order);
                        plot_counter = plot_counter + 1;
                        hold on
                    case 5
                        colour = colors{4};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(cell2mat(plot_structure(5, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("D:%d", Order);
                        plot_counter = plot_counter + 1;
                        hold on
                    otherwise
                        error("unknown type")
                end
            end

            title(simul.title)
            xlabel('$  \beta \cdot J$', 'Interpreter', 'latex', 'FontSize', 11)
            ylabel('$  \epsilon $', 'Interpreter', 'latex', 'FontSize', 11)
            legend(legend_Arr, 'Location', 'northwest', 'NumColumns', 2)

            ylim([0, 10])

            figure(gcf)

        end

        hold off

        filename = sprintf('../../auto_fig/comp%s.pdf', datestr(now, 'mm-dd-yy_HH-MM-SS'));
        saveas(gcf, filename)

    end
end
