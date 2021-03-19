function x_cell = solve_non_lin(obj, root_patterns, extended_patterns, pattern_root, pattern_permutations, maps, targets, con_cells, opts, lnprefact)
    p = inputParser;
    addParameter(p, 'optimise', 1)
    addParameter(p, 'Gradient', true)
    addParameter(p, 'maxit', 200)
    addParameter(p, 'Display', 'None')
    addParameter(p, 'PlotFcn', 'None')
    addParameter(p, 'Algoritm', 'levenberg-marquardt')
    %addParameter(p, 'Algoritm', 'trust-region')
    parse(p, opts)

    if nargin < 7
        lnprefact = obj.nf;
    end

    options = optimoptions('fsolve', 'Display', p.Results.Display, ...
        'Algorithm', p.Results.Algoritm, ...
        'MaxIterations', p.Results.maxit, ...
        'SpecifyObjectiveGradient', p.Results.Gradient, ...
        'FunctionTolerance', 1e-30, ...
        'StepTolerance', 1e-15, ...
        'PlotFcn', 'optimplotfirstorderopt', ...
        'OptimalityTolerance', 1e-40); %for trust region

    %

    num_patterns = size(root_patterns, 2);

    x_sizes = cell(1, num_patterns);
    begin_vec = [];

    for i = 1:size(root_patterns, 2)
        pattern = root_patterns{i};

        if size(pattern, 2) ~= 4
            error("fetch boundary matrices here")
        end

        tens = obj.PEPO_cell{pattern(1) + 1, pattern(2) + 1, pattern(3) + 1, pattern(4) + 1};
        x_sizes{i} = size(tens);
        begin_vec = [begin_vec, reshape(tens, 1, [])];
    end

    if p.Results.optimise == 1
        all_patterns = [root_patterns, extended_patterns];
        [con_cells, targets] = optimize_con_cells(obj, maps, con_cells, all_patterns, targets, lnprefact);
    end

    if p.Results.Algoritm == "trust-region"
        [~, g] = get_value_and_grad(obj, maps, con_cells, root_patterns, extended_patterns, pattern_root, pattern_permutations, targets, begin_vec, x_sizes, lnprefact);
        j_pat = (g ~= 0) * 1;
        options = optimoptions(options, 'JacobPattern', j_pat, 'SubproblemAlgorithm', 'cg');
    end

    x = fsolve(@(x) get_value_and_grad(obj, maps, con_cells, root_patterns, extended_patterns, pattern_root, pattern_permutations, targets, x, x_sizes, lnprefact), begin_vec, options);

    x_cell = split_x (x, x_sizes);

    for i = 1:numel(extended_patterns)
        real_pat = [1, 2, pattern_permutations{i} + 2];
        x_cell = [x_cell, permute(x_cell{pattern_root(i)}, real_pat)];
    end

    %             if obj.testing == 1
    %                 for i = 1:size(maps,2)
    %                     F = obj.get_value_and_grad(maps(i), con_cells2(i),root_patterns,targets2(i),begin_vec,x_sizes);
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
