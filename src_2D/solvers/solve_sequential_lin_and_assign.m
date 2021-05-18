function [obj, ln_prefact_out, err] = solve_sequential_lin_and_assign(obj, map, patterns, ln_prefact, opts, extended_patterns_permutations)
    %solves the given problem by repeatedly solving the different tensors appearing in the patterns using the linear solver.
    %assignfn -> obj = assignfn(obj, pattern) copies the given pattern to al equivalent patterns (e.g. rotated version)

    obj = fill_rand(obj, patterns, 1 / exp(obj.nf));

    if nargin >= 6
        assert(numel(patterns) == 1)

        [extended_patterns, ~, ~] = extend_pattern(patterns, extended_patterns_permutations);
        assignfn = @(x) assign_perm(x, patterns{1});

        all_patterns = [patterns, extended_patterns];
        obj = assignfn(obj);
    else
        all_patterns = patterns;
    end

    p = inputParser;
    addParameter(p, 'maxit', 200)
    addParameter(p, 'solved', 1e-14)
    addParameter(p, 'minchange', 1 - 1e-4)
    addParameter(p, 'display', 0)
    addParameter(p, 'counter', 20)
    addParameter(p, 'step', 1)
    parse(p, opts);

    counter = p.Results.counter;
    step = p.Results.step;

    maxit = p.Results.maxit;

    %precontract for speeed with matrix calculation and substract pattern
    %cells

    [all_con_cells, pat_cells] = get_valid_contractions(obj, map, struct('pattern', {patterns}));

    con_cells = all_con_cells(pat_cells);
    pat_cells = pat_cells(pat_cells == 1);

    assert(sum(pat_cells) == 1);

    obj = cell2matrix(obj);
    [target, ln_prefact_out] = H_exp(obj, map, ln_prefact, true);

    target_orig = target;

    target = target - contract_network(obj, map, struct('matrix', 1, 'lnprefact', ln_prefact_out));
    target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(obj.dim^2, map.N));

    target_site = -contract_con_cells(obj, map, ln_prefact_out, -target_site, con_cells);

    target_orig = target_orig - permute(reshape(target_site, dimension_vector(obj.dim, 2 * map.N)), [1:2:2 * map.N, 2:2:2 * map.N]);

    %find occuring patterns (needed for rescaling)
    num_pats = numel(all_patterns);
    nums = zeros(num_pats, 1) -1;

    for pat = 1:num_pats
        for i = 1:map.N
            if same_pattern(con_cells{1}{1}{i}, all_patterns{pat})
                nums(pat) = i;
                break;
            end
        end
    end

    rescale_patters = all_patterns(nums ~= -1);

    %
    prev_err = Inf; %calculate_error(obj, map.num_map, obj.numopts);
    for j = 1:maxit
        for i = 1:numel(patterns)
            if nargin >= 6
                leg = patterns{1};
                prev_A = obj.PEPO_cell{leg(1) + 1, leg(2) + 1, leg(3) + 1, leg(4) + 1};
            end

            [obj, ~, err] = solve_lin_and_assign(obj, map, patterns(i), ln_prefact_out, struct('target_site', target_site, 'all_con_cells', {con_cells}, 'pat_cells', pat_cells));
        end

        %rescale patterns such that they have same size
        obj = rescale_PEPO_pattern(obj, rescale_patters);

        if nargin >= 6

            leg = patterns{1};
            new = obj.PEPO_cell{leg(1) + 1, leg(2) + 1, leg(3) + 1, leg(4) + 1};
            old = prev_A;

            obj.PEPO_cell{leg(1) + 1, leg(2) + 1, leg(3) + 1, leg(4) + 1} = old + step * (new - old);

            obj = assignfn(obj);

            err_vect = contract_con_cells(obj, map, ln_prefact_out, target_site, con_cells);
            a = permute(reshape(err_vect, dimension_vector(obj.dim, 2 * map.N)), [1:2:2 * map.N, 2:2:2 * map.N]);

            err = calculate_error_core(a, target_orig);
        end

        if p.Results.display == 1
            fprintf("%3d: %.4e\n", j, err);
        end

        if err < p.Results.solved
            if p.Results.display == 1
                fprintf('under threshold')
            end
            break;
        end

        meas = err / prev_err;

        if err / prev_err > p.Results.minchange
            counter = counter - 1;

            if counter < 0
                if obj.testing == 1
                    fprintf("not solved, res err %.4e\n", err);
                end
                break;
            end

            if nargin >= 6
                if meas > 1 %step to big
                    step = step * 0.9;
                else
                    step = step * 1.1;
                end

                if obj.testing == 2
                    disp(step)
                end

            end
        end

        prev_err = err;
    end
end
