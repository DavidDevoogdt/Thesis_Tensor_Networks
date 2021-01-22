function x_cell = svd_x_cell(x,dims,bond_pairs,nums)
    %split x in different cells. SVD across bonds

    switch numel(bond_pairs)
        case 0
            x_cell = {reshape(x, dims{1})};
        case 1
            pair = bond_pairs{1};

            num1 = pair(1);
            num2 = pair(2);

            i1 = find(nums == num1);
            i2 = find(nums == num2);

            dims1 = dims{i1};
            dims2 = dims{i2};

            mask1 = dims1 == -1;
            mask2 = dims2 == -1;

            d1 = prod(dims1(~mask1));
            d2 = prod(dims2(~mask2));

            assert(d1 == d2)

            n1 = find(mask1);
            n2 = find(mask2);

            dim1_alt = [prod(dims1(1:n1 - 1)), prod(dims1(n1 + 1:end))];
            dim2_alt = [prod(dims2(1:n2 - 1)), prod(dims2(n2 + 1:end))];

            [U, S, V] = svd(reshape(x, d1, d1));

            sqrt_S = diag(diag(S).^0.5);

            l = permute(reshape(U * sqrt_S, dim1_alt(1), dim1_alt(2), []), [1, 3, 2]);
            r = permute(reshape(sqrt_S * V', [], dim2_alt(1), dim2_alt(2)), [2, 1, 3]);

            dims1(mask1) = size(l, 2);
            dims2(mask2) = size(r, 2);

            x_cell{i1} = reshape(l, dims1);
            x_cell{i2} = reshape(r, dims2);
        otherwise
            error("not implemted")
    end
end