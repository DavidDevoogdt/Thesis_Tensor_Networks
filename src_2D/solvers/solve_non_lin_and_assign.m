function [obj, ln_prefact] = solve_non_lin_and_assign(obj, maps, root_patterns, ln_prefact, opts, extended_patterns_permutations)

    d = obj.dim;

    if nargin < 6
        extended_patterns_permutations = cell(size(root_patterns));
    end

    [extended_patterns, pattern_root, pattern_permutations] = extend_pattern(root_patterns, extended_patterns_permutations);

    obj = fill_rand(obj, root_patterns);
    obj = fill_rand(obj, extended_patterns);

    all_patterns = [root_patterns, extended_patterns];

    nmaps = numel(maps);
    targets = cell(1, nmaps);
    con_cells = cell(1, nmaps);

    ln_prefact_out = -1;

    for i = 1:nmaps
        if ln_prefact_out == -1
            [target, ln_prefact_out] = H_exp(obj, maps{i}, ln_prefact, true);
        else
            [target, ~] = H_exp(obj, maps{i}, ln_prefact_out, false);
        end
        targets{i} = reshape(permute(target, site_ordering_permute(maps{i}.N)), dimension_vector(d^2, maps{i}.N));
        con_cells{i} = get_valid_contractions(obj, maps{i}, struct('max_index', obj.current_max_index, 'pattern', {all_patterns}));
    end

    mul_factor = exp(ln_prefact_out - obj.nf);

    %refill with new prefact
    %fprintf("old nf %.4e new %.4e",ln_prefact,ln_prefact_out);

    init_val = 1e-3 / mul_factor;

    obj = fill_rand(obj, root_patterns, init_val, true);
    obj = fill_rand(obj, extended_patterns, init_val, true);

    %     if nmaps ==1
    %        maps = {maps};
    %     end

    x_cell = solve_non_lin(obj, root_patterns, extended_patterns, pattern_root, pattern_permutations, maps, targets, con_cells, opts, ln_prefact_out);

    for i = 1:size(all_patterns, 2)
        obj.PEPO_cell{all_patterns{i}(1) + 1, all_patterns{i}(2) + 1, all_patterns{i}(3) + 1, all_patterns{i}(4) + 1} = x_cell{i} * (mul_factor);
        %fprintf("%.4e ", max(abs(reshape(x_cell{i} * mul_factor, [], 1))));
        if obj.testing == 2
            fprintf("%.4e ", norm(reshape(x_cell{i} * mul_factor, [], 1), 2));
        end
    end

end
