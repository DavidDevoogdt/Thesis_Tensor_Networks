function [A, x0_shape] = differentiate_x(obj, num, map, con_cell, ln_prefactor, x, pattern, extended_patterns, pattern_root, pattern_permutations)
    % takes map with a number of cells removed, and calculates the the dervative
    % if only 1 element is removed, the output is reordered to calculate the gradient easily
    % for linear problem, calculates A x = B
    % x, pattern, extended_patterns, pattern_root, pattern_permutations all represent the for current point. get_value and grad uses this
    % to calculate value at a given point

    if nargin < 5
        ln_prefactor = obj.nf;
    end

    seq = map.seq;
    final_order = map.final_order;
    leg_list = map.leg_list;

    mask = map.leg_list_mask;
    num_removed = sum(~mask);

    legs = con_cell{1};

    if nargin < 6
        temp_list = fetch_PEPO_cells(obj, map, legs, ln_prefactor);
    else
        temp_list = fetch_PEPO_cells(obj, map, legs, ln_prefactor, pattern, x, extended_patterns, pattern_root, pattern_permutations);
    end

    x0 = temp_list(~mask);
    x0_shape = size(x0);

    temp_list(~mask) = []; %remove x0
    leg_list = leg_list(mask);

    if map.N - num_removed == 0
        A = 1;
    else
        A = ncon(temp_list, leg_list, seq, final_order);
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
            %put physical dims per site, leave external legs alone
            perm_vector = [site_ordering_permute(map.N2 - num_removed); ((2 * (map.N2 - num_removed) + 1):size(size(A), 2)).'];

            A = permute(A, perm_vector);
        end
    end
end
