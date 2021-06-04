function test

    %compare_models(["t_ising"], 1i*10.^(-3:0.5:1.5), [2, 3, 4], [3, 4, 5,
    %6],1) %imaginary
    
    compare_models(["t_ising"], 10.^(-3:0.2:1.5), [2], [3,5,7],0)
    
    %compare_M( ["t_ising","Heisenberg_2D"] , 10.^(-3:0.9:1),  [2,3,4,5,6,7,8,9] );

end

function compare_models(simulatiemodellen, beta_arr, types, order,time)
    models_len = size(simulatiemodellen, 2);

    for round = 1:models_len
        %change this
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        model = models(simulatiemodellen(round), struct);

        simul.Order_arr = order;
        simul.types = types;

        simul.M = 10;
        simul.beta_arr = beta_arr;
        simul.cyclic = 1;

        map = create_map(1:simul.M, struct("numbered", true, "h_cyclic", simul.cyclic));
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
        
        small = 1;
        
        if small == 0
            ncols =2;
            x_width = 15;
            y_width = 10;
        else
            ncols =1;
            x_width = 10;
            y_width = 8;
        end
        
        set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, x_width, y_width], 'PaperSize', [x_width, y_width])
        %%%

        opts = struct('max_bond_dim',40,'complex', false);
        %opts = struct('complex', false);
        opts.inv_eps = 1e-12;
        %opts = struct;

        %loop over different orders
        for j = 1:order_size
            opts.order = simul.Order_arr(j);
            plot_structure = cell(5, beta_len);

            %loop over temps
            for i = 1:beta_len
                opts.beta = simul.beta_arr(i);

                 if time == 1
                     fprintf("M %d time %.4e order %d", simul.M, imag(opts.beta), opts.order)
                 else
                     
                    fprintf("M %d beta %.4e order %d", simul.M, opts.beta, opts.order);
                 end

                %loop of simulation types
                for t = 1:size(simul.types, 2)
                    switch simul.types(t)
                        case 2
                            handle = @make_PEPO_1D_type_A;
                            pepo = PEPO(model, opts, handle);

                            fprintf(".");

                            err_02 = calculate_error(pepo, map, [], 1);

                            fprintf(" err 02 %.4e", err_02);
                            plot_structure{2, i} = err_02;
                        case 3
                            handle = @make_PEPO_1D_type_E;
                            pepo = PEPO(model, opts, handle);

                            fprintf(".");
                            err_03 = calculate_error(pepo, map, [], 1);

                            fprintf(" err 03 %.4e", err_03);
                            plot_structure{3, i} = err_03;
                        case 4
                            handle = @make_PEPO_1D_type_F;
                            pepo = PEPO(model, opts, handle);
                            fprintf(".");
                            err_04 = calculate_error(pepo, map, [], 1);

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

            x_axis = abs(simul.beta_arr);
            
            %plotting loop for current order
            for t = 1:size(simul.types, 2)
                switch simul.types(t)
                    case 2
                        colour = colors{1};
                        colour(4) = alphas(j);
                        loglog(x_axis, abs(cell2mat(plot_structure(2, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("A:%d", opts.order);
                        plot_counter = plot_counter + 1;
                        hold on
                    case 3
                        colour = colors{2};
                        colour(4) = alphas(j);
                        loglog(x_axis, abs(cell2mat(plot_structure(3, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("E:%d", opts.order);
                        plot_counter = plot_counter + 1;
                        hold on
                    case 4
                        colour = colors{3};
                        colour(4) = alphas(j);
                        loglog(x_axis, abs(cell2mat(plot_structure(4, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("F:%d", opts.order);
                        plot_counter = plot_counter + 1;
                        hold on

                    case 5
                        colour = colors{4};
                        colour(4) = alphas(j);
                        loglog(x_axis, abs(cell2mat(plot_structure(5, :))), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("D:%d", opts.order);
                        plot_counter = plot_counter + 1;
                        hold on
                    otherwise
                        error("unknown type")
                end
            end

            title(model.title)
            if time == 1
                xlabel('$  t \cdot J$', 'Interpreter', 'latex', 'FontSize', 11)
            else
                xlabel('$  \beta \cdot J$', 'Interpreter', 'latex', 'FontSize', 11)
            end
            
            ylabel('$  \epsilon $', 'Interpreter', 'latex', 'FontSize', 11)
            legend(legend_Arr, 'Location', 'northwest', 'NumColumns', ncols)

            ylim([0, 10])

            figure(gcf)

        end

        hold off

        filename = sprintf('./auto_fig/comp%s.pdf', datestr(now, 'mm-dd-yy_HH-MM-SS'));
        saveas(gcf, filename)

    end
end

function compare_M(simulatiemodellen, beta_arr, M)
    models_len = size(simulatiemodellen, 2);

    for round = 1:models_len
        %change this
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        model = models(simulatiemodellen(round), struct);

        opts.order = 5;
        opts.testing = 1;

        simul.types = [2, 3];

        simul.M = M;
        simul.beta_arr = beta_arr;
        simul.cyclic = 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%plot stuff
        M_size = size(simul.M, 2);
        line_spec = ["-", "--", "-.", ":", "-", "--", "-.", ":"];
        alphas = [1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5];
        colors = {[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 1]}; %reserved for specific type

        legend_Arr = cell(size(simul.types, 2) * M_size, 1);
        legend_Arr(:) = {"todo"};
        beta_len = size(simul.beta_arr, 2);

        plot_counter = 1;
        %hold off
        figure();
        x_width = 15;
        y_width = 10;
        set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, x_width, y_width], 'PaperSize', [x_width, y_width])
        %%%

        %loop over different orders

        plot_structure = zeros(5, M_size, beta_len);

        for i = 1:beta_len
            opts.beta = simul.beta_arr(i);

            fprintf("beta %.4e order %d \n", opts.beta, opts.order);

            for t = 1:size(simul.types, 2)
                switch simul.types(t)
                    case 2
                        handle = @make_PEPO_1D;
                        pepo02 = PEPO(model, opts, handle);
                    case 3
                        handle = @make_PEPO_1D_type_E;
                        pepo03 = PEPO(model, opts, handle);

                    case 4
                        handle = @make_PEPO_1D_G;
                        pepo04 = PEPO(model, opts, handle);
                    case 5
                        error("")
                    otherwise
                        error("unknown type")
                end
            end

            for j = 1:M_size

                map = create_map(1:M(j), struct("numbered", true, "h_cyclic", simul.cyclic));

                for t = 1:size(simul.types, 2)
                    switch simul.types(t)
                        case 2
                            err_02 = calculate_error(pepo02, map, [], 1);
                            %fprintf(" err 02 %.4e", err_02);
                            plot_structure(2, j, i) = err_02;
                        case 3
                            err_03 = calculate_error(pepo03, map, [], 1);

                            %fprintf(" err 03 %.4e", err_03);
                            plot_structure(3, j, i) = err_03;
                        case 4
                            err_04 = calculate_error(pepo04, map, [], 1);
                            %fprintf(" err 04 %.4e", err_04);
                            plot_structure(4, j, i) = err_04;

                        case 5
                            error("")
                        otherwise
                            error("unknown type")
                    end
                end
            end

        end

        for j = 1:M_size

            for t = 1:size(simul.types, 2)
                switch simul.types(t)
                    case 2
                        colour = colors{1};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(reshape(plot_structure(2, j, :), [], 1)), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("A: M %d", M(j));
                        plot_counter = plot_counter + 1;
                        hold on
                    case 3
                        colour = colors{2};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(reshape(plot_structure(3, j, :), [], 1)), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("E: M %d", M(j));
                        plot_counter = plot_counter + 1;
                        hold on
                    case 4
                        error('check this')
                        colour = colors{3};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(reshape(plot_structure(4, j, :), [], 1)), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("A: M %d", M(j));
                        plot_counter = plot_counter + 1;
                        hold on
                    case 5
                        error('check this')
                        colour = colors{4};
                        colour(4) = alphas(j);
                        loglog(simul.beta_arr, abs(reshape(plot_structure(5, j, :), [], 1)), "LineStyle", line_spec(j), "Color", colour);
                        legend_Arr{plot_counter} = sprintf("D: M %d", M(j));
                        plot_counter = plot_counter + 1;
                        hold on
                    otherwise
                        error("unknown type")
                end
            end

        end

        title(model.title)
        xlabel('$  \beta \cdot J$', 'Interpreter', 'latex')
        ylabel('$  \epsilon $', 'Interpreter', 'latex')
        legend(legend_Arr, 'Location', 'northwest', 'NumColumns', 2)

        ylim([0, 10])

        figure(gcf)

        hold off

        filename = sprintf('./auto_fig/comp%s.pdf', datestr(now, 'mm-dd-yy_HH-MM-SS'));
        saveas(gcf, filename)

    end
end
