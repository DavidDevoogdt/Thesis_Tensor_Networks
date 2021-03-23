function x_cell = svd_x_cell(x, dims, bond_pairs, nums, split_dim)
    if nargin < 5
        split_dim = -1;
    end

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

            %assert(d1 == d2)

            n1 = find(mask1);
            n2 = find(mask2);

            dim1_alt = [prod(dims1(1:n1 - 1)), prod(dims1(n1 + 1:end))];
            dim2_alt = [prod(dims2(1:n2 - 1)), prod(dims2(n2 + 1:end))];

            x_res = reshape(x, d1, d2);

            [U, S, V] = svd(x_res);

            %err = U*S*V'-x_res;
            %err = U*S_l*S_r*V'-x_res;

            if split_dim == -1

                if d1 == d2
                    split_dim = d1;
                else
                    error("provide split dimension")
                end
            end

            PU = eye(size(U, 1), split_dim);
            PR = eye(split_dim, size(V, 1));
            DS = diag(S);
            sqrt_S = diag(DS(1:split_dim).^0.5);

            L = U * PU * sqrt_S;
            R = sqrt_S * PR * V';

            err = L * R - x_res;

            parity = find(mask1) > 4; %order of multiplication

            if parity == 1 %in right order
                l = permute(reshape(L, dim1_alt(1), dim1_alt(2), []), [1, 3, 2]);
                r = permute(reshape(R, [], dim2_alt(1), dim2_alt(2)), [2, 1, 3]);

                dims1(mask1) = split_dim;
                dims2(mask2) = split_dim;

                x_cell{i1} = reshape(l, dims1);
                x_cell{i2} = reshape(r, dims2);

            else
                l = permute(reshape(L, dim2_alt(1), dim2_alt(2), []), [1, 3, 2]);
                r = permute(reshape(R, [], dim1_alt(1), dim1_alt(2)), [2, 1, 3]);

                dims1(mask1) = split_dim;
                dims2(mask2) = split_dim;

                x_cell{i1} = reshape(r, dims1);
                x_cell{i2} = reshape(l, dims2);

            end

        otherwise
            error("not implemented")
    end
end
