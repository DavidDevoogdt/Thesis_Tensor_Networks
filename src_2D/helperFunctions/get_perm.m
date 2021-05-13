function dd = get_perm(pat, mask, rot)
    %rot mask: 0 -> all perm
    % 1-> keep rotation order

    if nargin < 2
        mask = [0, 0, 0, 0];
    end

    orig_perms = 1:4;

    nr = mask ~= 1;

    all_perms = perms(orig_perms(nr));
    perm_pat = pat(all_perms);
    [~, ia, ~] = unique(perm_pat, 'rows');

    perm_vect = repmat(orig_perms, numel(ia), 1);
    perm_vect(:, nr) = all_perms(ia, :);

    perm_val = repmat(pat, numel(ia), 1);
    perm_val(:, nr) = perm_pat(ia, :);

    if sum(mask == 1) ~= 0 %also include  rotations
        pvect = perm_vect;
        pval = perm_val;

        for k = 1:3
            perm_val = cat(1, perm_val, circshift(pval, k, 2));
            perm_vect = cat(1, perm_vect, circshift(pvect, k, 2));
        end
    end

    %filter in groyps according to left and right matrix inverse

    A = unique(sort(perm_val(:, 1:2)')', 'rows');
    B = unique(sort(perm_val(:, 3:4)')', 'rows');

    perm_cell = cell(size(A, 1), size(B, 1), 2);

    for i = 1:size(perm_val, 1)
        a = find(ismember(A, sort(perm_val(i, 1:2)), 'rows'));
        b = find(ismember(A, sort(perm_val(i, 3:4)), 'rows'));

        perm_cell{a, b, 1} = [perm_cell{a, b, 1}; perm_val(i, :)];
        perm_cell{a, b, 2} = [perm_cell{a, b, 2}; perm_vect(i, :)];
    end

    nempty = cellfun(@(x)~isempty(x), perm_cell);

    dd = reshape(perm_cell(nempty), [], 2);
end
