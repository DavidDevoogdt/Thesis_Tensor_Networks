function [x_cell, con_cells, targets] = solve_non_lin(obj, root_patterns, extended_patterns, pattern_root, pattern_permutations, maps, targets, con_cells, opts, lnprefact)
    %do not use directly, but through solve_non_lin_and_assign
    %opts: see below

    p = inputParser;
    addParameter(p, 'optimise', 1)
    addParameter(p, 'Gradient', true)
    addParameter(p, 'maxit', 400)
    addParameter(p, 'Display', 'None')
    addParameter(p, 'PlotFcn', []) %'optimplotfirstorderopt'
    addParameter(p, 'Algoritm', 'levenberg-marquardt') %"trust-region"
    %addParameter(p, 'Algoritm', 'trust-region')
    parse(p, opts)

    if nargin < 7
        lnprefact = obj.nf;
    end

    options = optimoptions('fsolve', 'Display', p.Results.Display, ...
        'Algorithm', p.Results.Algoritm, ...
        'MaxIterations', p.Results.maxit, ...
        'SpecifyObjectiveGradient', p.Results.Gradient, ...
        'FunctionTolerance', 1e-15, ...
        'StepTolerance', 1e-13, ...
        'PlotFcn', p.Results.PlotFcn, ...
        'InitDamping', 1e-5, ...
        'OptimalityTolerance', 1e-20); %for trust region

    num_patterns = size(root_patterns, 2);

    x_sizes = cell(1, num_patterns);
    begin_vec = [];

    %assemble beginvalue for fsovle
    for i = 1:size(root_patterns, 2)
        pattern = root_patterns{i};

        if size(pattern, 2) ~= 4
            error("fetch boundary matrices here")
        end

        tens = obj.PEPO_cell{pattern(1) + 1, pattern(2) + 1, pattern(3) + 1, pattern(4) + 1};
        x_sizes{i} = size(tens);
        begin_vec = [begin_vec, reshape(tens, 1, [])];
    end

    %presimplify contraction each iteration
    if p.Results.optimise == 1
        all_patterns = [root_patterns, extended_patterns];
        [con_cells, targets] = optimize_con_cells(obj, maps, con_cells, all_patterns, targets, lnprefact);
    end

    %pre compute jacobian pattern
    if p.Results.Algoritm == "trust-region"
        [~, g] = get_value_and_grad(obj, maps, con_cells, root_patterns, extended_patterns, pattern_root, pattern_permutations, targets, begin_vec, x_sizes, lnprefact);
        j_pat = (g ~= 0) * 1;
        options = optimoptions(options, 'JacobPattern', j_pat, 'SubproblemAlgorithm', 'cg');
    end

    x = fsolve(@(x) get_value_and_grad(obj, maps, con_cells, root_patterns, extended_patterns, pattern_root, pattern_permutations, targets, x, x_sizes, lnprefact), begin_vec, options);

    %assign output to tensors
    x_cell = split_x (x, x_sizes);

    for i = 1:numel(extended_patterns)
        real_pat = [1, 2, pattern_permutations{i} + 2];
        x_cell = [x_cell; permute(x_cell{pattern_root(i)}, real_pat)];
    end

end
