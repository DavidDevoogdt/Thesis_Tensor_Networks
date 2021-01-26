function x_cell = solve_non_lin(obj, patterns, maps, targets, con_cells, opts, lnprefact)
    p = inputParser;
    addParameter(p, 'optimise', 1)
    addParameter(p, 'Gradient', true)
    parse(p, opts)

    if nargin < 7
        lnprefact = obj.nf;
    end

    options = optimoptions('fsolve', 'Display', 'iter-detailed', ...
        'Algorithm', 'levenberg-marquardt', ...
        ... %'MaxIterations',2000,...
        'SpecifyObjectiveGradient', p.Results.Gradient, ...
        'FunctionTolerance', 1e-40, ...
        'StepTolerance', 1e-10);

    %'PlotFcn', 'optimplotfirstorderopt',...

        num_patterns = size(patterns, 2);

    x_sizes = cell(1, num_patterns);
    begin_vec = [];

    for i = 1:size(patterns, 2)
        pattern = patterns{i};

        if size(pattern, 2) ~= 4
            error("fetch boundary matrices here")
        end

        tens = obj.PEPO_cell{pattern(1) + 1, pattern(2) + 1, pattern(3) + 1, pattern(4) + 1};
        x_sizes{i} = size(tens);
        begin_vec = [begin_vec, reshape(tens, 1, [])];
    end

    if p.Results.optimise == 1
        [con_cells, targets] = optimize_con_cells(obj, maps, con_cells, patterns, targets, lnprefact);
    end

    x = fsolve(@(x) get_value_and_grad(obj, maps, con_cells, patterns, targets, x, x_sizes, lnprefact), begin_vec, options);

    x_cell = split_x (x, x_sizes);

    %             if obj.testing == 1
    %                 for i = 1:size(maps,2)
    %                     F = obj.get_value_and_grad(maps(i), con_cells2(i),patterns,targets2(i),begin_vec,x_sizes);
    %                     norm = sum(F.^2)^0.5
    %                 end
    %
    %                 for i = 1:size(tf,1)
    %                     A = 0;
    %                     for j = 1:size(tf{i},2)
    %                         tensors = obj.fetch_PEPO_cells(map_2, tf{i}{j}{1} );
    %                         vect  =  ncon( tensors, map_2.leg_list );
    %                         leg_1_size =  prod(size(vect,5:8));
    %                         vect = reshape(vect,d^2,d^2, leg_1_size,[]  );
    %                         perm = reshape( permute(vect, [3,1,2,4]), leg_1_size*d^2,[]);
    %                         A=A+perm;
    %
    %                     end
    %                     S= (svds(A,2).^0.5)
    %                 end
    %
    %             end

end
