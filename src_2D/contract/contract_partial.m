function [A, x0_shape] = contract_partial(obj, num, map, con_cells, ln_prefactor, x, pattern)

    %patterns: match x with pattern given pattern for subs
    if nargin < 5
        ln_prefactor = obj.nf;
    end

    seq = map.seq;
    final_order = map.final_order;
    leg_list = map.leg_list;

    % do actual contractions
    A = [];

    mask = map.leg_list_mask;
    num_removed = sum(~mask);

    for con_cell_index = 1:size(con_cells, 2)
        legs = con_cells{con_cell_index}{1};

        if nargin < 6
            temp_list = fetch_PEPO_cells(obj, map, legs, ln_prefactor);
        else
            temp_list = fetch_PEPO_cells(obj, map, legs, ln_prefactor, pattern, x);
        end

        x0 = temp_list(~mask);
        x0_shape = size(x0);

        temp_list(~mask) = []; %remove x0
        leg_list = leg_list(mask);

        if map.N - num_removed == 0
            T = 1;
        else
            T = ncon(temp_list, leg_list, seq, final_order);
        end

        if isempty(A)
            A = T;
        else
            A = A + T;
        end

    end

    if num_removed == 1

        %reorder such that gradient is easy to compute
        %format: [physical indices before,1,physical indices after, external legs before,1,external legs after, bonds to x ]
        if map.is_x_border(num) || map.is_x_border(num)
            perm_vector = [site_ordering_permute(map.N2); ((2 * map.N2 + 1):size(size(A), 2)).'];
            A = reshape(permute(A, perm_vector), [], prod(x0_shape(1:2)));
        else
            a_size = size(A);

            x0_list = map.leg_list{num};

            ext_legs = x0_list(3:end);
            smallest_ind = min(-ext_legs(ext_legs < 0));

            num_removed = sum(~map.leg_list_mask);
            num1 = 2 * (map.N2 - num_removed);

            if ~isempty(smallest_ind)
                idx = find(final_order == -(smallest_ind -1));
            else
                idx = num1;
            end

            num2 = idx;
            num3 = size(size(A), 2) - map.ii;

            size1 = a_size(1:num1); %ij indices
            size2 = a_size(num1 + 1:num2); %external legs before
            size3 = a_size(num2 + 1:num3); %external legs after
            size4 = a_size(num3 + 1:end); %bond to x0

            perm_vector = [site_ordering_permute(map.N2 - 1); ((2 * map.N2 - 1):size(size(A), 2)).'];
            A = reshape(permute(A, perm_vector), prod(size1), prod(size2), 1, prod(size3), prod(size4));
        end

    else

        if map.N - num_removed ~= 0
            perm_vector = [site_ordering_permute(map.N2 - num_removed); ((2 * (map.N2 - num_removed) + 1):size(size(A), 2)).'];

            A = permute(A, perm_vector);
        end
    end
end
