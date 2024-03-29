function obj = assign_perm(obj, pat, rot_mask)

    if nargin < 3
        rot_mask = [0, 0, 0, 0];
    end

    orig_perms = 1:4;

    nr = rot_mask ~= 1;

    %orig = [a, a, 0, 0];
    all_perms = perms(orig_perms(nr));
    perm_pat = pat(all_perms);
    [~, ia, ~] = unique(perm_pat, 'rows');

    perm_vect = repmat(orig_perms, numel(ia), 1);
    perm_vect(:, nr) = all_perms(ia, :);

    perm_val = repmat(pat, numel(ia), 1);
    perm_val(:, nr) = perm_pat(ia, :);

    if sum(rot_mask) ~= 0 %also include  rotations
        pvect = perm_vect;
        pval = perm_val;

        for k = 1:3
            perm_val = cat(1, perm_val, circshift(pval, k, 2));
            perm_vect = cat(1, perm_vect, circshift(pvect, k, 2));
        end
    end

    pp = pat + [rot_mask == 2] .* [0, 0, 1, 1];

    tensor0 = obj.PEPO_cell{pp(1) + 1, pp(2) + 1, pp(3) + 1, pp(4) + 1};

    tensor = sym_com(tensor0, pat, rot_mask);

    for i = 1:size(perm_vect, 1)
        pv = perm_vect(i, :);
        rp = rot_mask(perm_vect(i, :));

        vect = [1, 2, pv + 2];
        rpat = perm_val(i, :);

        rpat2 = rpat + [rp == 2] .* [0, 0, 1, 1];
        rpat2 = rpat2 + [rp == -2] .* [1, 1, 0, 0];

        eq_tensor = permute(tensor, vect);

        %if isempty(obj.PEPO_cell{rpat(1) + 1, rpat(2) + 1, rpat(3) + 1, rpat(4) + 1})
        obj.PEPO_cell{rpat(1) + 1, rpat(2) + 1, rpat(3) + 1, rpat(4) + 1} = [];
        obj.PEPO_cell{rpat2(1) + 1, rpat2(2) + 1, rpat2(3) + 1, rpat2(4) + 1} = eq_tensor;
        %end

    end
end

function tens = sym_com(tens, pat, rot_mask)
    unique_nums = unique(pat(~rot_mask));
    orig_perm = 1:4;

    orig_tens = tens;

    for i = 1:numel(unique_nums)
        n = unique_nums(i);
        %elses(i);

        mask = find(pat == n);

        tens0 = tens;
        tens = zeros(size(tens));

        [all_perms] = perms(mask);
        np = size(all_perms, 1);

        for ii = 1:np
            this_perm = orig_perm;
            this_perm(mask) = all_perms(ii, :);

            this_perm = [1, 2, this_perm + 2];

            tens = tens + permute(tens0, this_perm) ./ np;
        end

    end

    %sym_err =svds( reshape( orig_tens-tens, [],1)  )
end
