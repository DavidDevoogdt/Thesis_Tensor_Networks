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

            if d1 ~= d2
                if split_dim ~= -1
                    error("not implemented")
                end

                if d1 > d2
                    L = U * S;
                    R = V';
                else
                    L = U;
                    R = S * V';
                end

            else
                sqrt_S = diag(diag(S).^0.5);
                s_dim = size(sqrt_S, 1);

                if split_dim ~= -1
                    sqrt_S_red = sqrt_S(1:split_dim, 1:split_dim);
                    S_l = eye(s_dim, split_dim) * sqrt_S_red;
                    S_r = sqrt_S_red * eye(split_dim, s_dim);

                    ds = diag(S);
                    err = sum(ds(split_dim + 1:end));

                else
                    S_l = sqrt_S;
                    S_r = sqrt_S;
                end
                L = U * S_l;
                R = S_r * V';
            end

            parity = mod(find(mask1), 2); %order of multiplication

            if parity == 1
                l = permute(reshape(L, dim1_alt(1), dim1_alt(2), []), [1, 3, 2]);
                r = permute(reshape(R, [], dim2_alt(1), dim2_alt(2)), [2, 1, 3]);

                dims1(mask1) = size(l, 2);
                dims2(mask2) = size(r, 2);

                x_cell{i1} = reshape(l, dims1);
                x_cell{i2} = reshape(r, dims2);

            else
                l = permute(reshape(L, dim2_alt(1), dim2_alt(2), []), [1, 3, 2]);
                r = permute(reshape(R, [], dim1_alt(1), dim1_alt(2)), [2, 1, 3]);

                dims1(mask1) = size(r, 2);
                dims2(mask2) = size(l, 2);

                x_cell{i1} = reshape(r, dims1);
                x_cell{i2} = reshape(l, dims2);

            end

        otherwise
            error("not implemented")
    end
end
