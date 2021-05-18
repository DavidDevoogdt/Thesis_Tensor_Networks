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
        model = models(simulatiemodellen(round), model_opts);

        simul.Order_arr = [5];
        simul.types = [2, 3];
        simul.M = 9;
        simul.beta_arr = 10.^(-1:0.1:0.5);
        %simul.beta_arr = 10.^(  -0:0.05:1);
        simul.cyclic = 1;

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

        opts = struct;

        %loop over different orders
        for j = 1:order_size
            opts.order = simul.Order_arr(j);
            plot_structure = cell(5, beta_len);

            %loop over temps
            for i = 1:beta_len
                opts.beta = simul.beta_arr(i);

                %mpo_base_matrix = mpo_base.H_exp(simul.M-1, 1, simul.cyclic); %buffered, not calculted again every time

                fprintf("M %d beta %.4e order %d", simul.M, opts.beta, opts.order);

                %loop of simulation types
                for t = 1:size(simul.types, 2)
                    switch simul.types(t)
                        case 2
                            handle = @make_PEPO_1D;
                            pepo = PEPO(model, opts, handle);

                            fprintf(".");
                            err_02 = calculate_error(pepo, 1:10, struct("numbered", true, "h_cyclic", simul.cyclic), 1);

                            fprintf(" err 02 %.4e", err_02);
                            plot_structure{2, i} = err_02;
                        case 3
                            handle = @make_PEPO_1D_double;
                            pepo = PEPO(model, opts, handle);

                            fprintf(".");
                            err_03 = calculate_error(pepo, 1:10, struct("numbered", true, "h_cyclic", simul.cyclic), 1);

                            fprintf(" err 03 %.4e", err_03);
                            plot_structure{3, i} = err_03;
                        case 4
                            handle = @make_PEPO_1D_G;
                            pepo = PEPO(model, opts, handle);
                            fprintf(".");
                            err_04 = calculate_error(pepo, 1:10, struct("numbered", true, "h_cyclic", simul.cyclic), 1);

                            fprintf(" err 04 %.4e", err_04);
                            plot_structure{4, i} = err_04;
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
                        legend_Arr{plot_counter} = sprintf("A:%d", opts.order);
                        plot_counter = plot_counter + 1;
                        hold on
                    case 4
                        colour = colors{3};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(cell2mat(plot_structure(4, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("A:%d", opts.order);
                        plot_counter = plot_counter + 1;
                        hold on
                    case 3
                        colour = colors{2};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(cell2mat(plot_structure(3, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("C:%d", opts.order);
                        plot_counter = plot_counter + 1;
                        hold on
                    case 5
                        colour = colors{4};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(cell2mat(plot_structure(5, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("D:%d", opts.order);
                        plot_counter = plot_counter + 1;
                        hold on
                    otherwise
                        error("unknown type")
                end
            end

            title(model.title)
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
