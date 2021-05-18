function [extened_patterns, pattern_root, pattern_permutations] = extend_pattern(patterns_in, extended_patterns_permutations)
    %get all patterns from given root pattern and its permutations

    n = numel(patterns_in);
    extened_patterns = {};
    pattern_root = [];

    pattern_permutations = {};

    for i = 1:n
        this_pat = patterns_in{i};
        perm = extended_patterns_permutations{i};
        m = numel(perm);
        for j = 1:m
            extened_patterns = [extened_patterns, this_pat(perm{j})];
            pattern_root = [pattern_root, i];
        end
        pattern_permutations = [pattern_permutations, perm];
    end
end
