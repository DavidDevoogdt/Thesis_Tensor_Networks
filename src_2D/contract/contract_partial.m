function [A, x0_shape] = contract_partial(obj, num, map, con_cells, ln_prefactor, x, pattern)

    %patterns: match x with pattern given pattern for subs

    if nargin < 5
        ln_prefactor = obj.nf;
    end

    seq = map.seq;
    final_order = map.final_order;
    leg_list = map.leg_list;

    % %
    % %             % create new contraction list without num and legs connected to
    % %             % num as last indices
    % %             l=map.h_bond_l_lookup{num};
    % %             r=map.h_bond_r_lookup{num};
    % %             u=map.v_bond_u_lookup{num};
    % %             d=map.v_bond_d_lookup{num};
    % %
    % %
    % %             con_list_cpy = map.leg_list;
    % %
    % %             ii = 0;
    % %
    % %             if ~isempty(l)
    % %                 ii = ii+1;
    % %                 pair = map.h_bonds{l};
    % %                 other = pair(1);
    % %
    % %                 if map.is_x_border(other)
    % %                     con_list_cpy{other}(2) = -(map.external_legs+ii);
    % %                 else
    % %                     con_list_cpy{other}(2+3) = -(map.external_legs+ii);
    % %                 end
    % %
    % %             end
    % %
    % %             if ~isempty(u)
    % %                 ii = ii+1;
    % %                 pair = map.v_bonds{u};
    % %                 other = pair(1);
    % %                 if map.is_y_border(other)
    % %                     con_list_cpy{other}(2) = -(map.external_legs+ii);
    % %                 else
    % %                     con_list_cpy{other}(2+4) = -(map.external_legs+ii);
    % %                 end
    % %             end
    % %
    % %             if ~isempty(r)
    % %                 ii = ii+1;
    % %                 pair = map.h_bonds{r};
    % %                 other = pair(2);
    % %                 if map.is_x_border(other)
    % %                     con_list_cpy{other}(1) = -(map.external_legs+ii);
    % %                 else
    % %                     con_list_cpy{other}(2+1) = -(map.external_legs+ii);
    % %                 end
    % %             end
    % %
    % %             if ~isempty(d)
    % %                 ii = ii+1;
    % %                 pair = map.v_bonds{d};
    % %                 other = pair(2);
    % %                 if map.is_y_border(other)
    % %                     con_list_cpy{other}(1) = -(map.external_legs+ii);
    % %                 else
    % %                     con_list_cpy{other}(2+2) = -(map.external_legs+ii);
    % %                 end
    % %             end
    % %
    % %             num_legs_cpy = con_list_cpy{num};
    % %
    % %             x0_list = con_list_cpy{ num};
    % %             con_list_cpy( num) = [];
    % %
    % %
    % %             final_order = -1:-1: -(map.external_legs+ii);
    % %             seq = 1:map.internal_legs; %contraction sequence
    % %
    % %
    % %             final_external_missing =sort( - num_legs_cpy(num_legs_cpy<0),'descend');
    % %             final_internal_missing = sort( num_legs_cpy(num_legs_cpy>0),'descend');
    % %
    % %             for l = 1:size(final_external_missing,2)
    % %                 final_order(  final_external_missing(l)  )= [];
    % %             end
    % %
    % %             for l = 1:size(final_internal_missing,2)
    % %                 seq(  final_internal_missing(l)  )= [];
    % %             end
    % %

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
            T = [1];
        else
            T = ncon(temp_list, leg_list, seq, final_order);
        end

        if isempty(A)
            A = T;
        else
            A = A + T;
        end

    end

    switch num_removed
        case 1
            %put in site ordering
            if map.is_x_border(num) || map.is_x_border(num)
                perm_vector = [site_ordering_permute(map.N2); ((2 * map.N2 + 1):size(size(A), 2)).'];
                A = reshape(permute(A, perm_vector), [], prod(x0_shape(1:2)));
            else
                a_size = size(A);

                x0_list = map.leg_list{num};

                ext_legs = x0_list(3:end);
                smallest_ind = min(-ext_legs(ext_legs < 0));

                idx = find(final_order == -(smallest_ind -1));

                num_removed = sum(~map.leg_list_mask);

                num1 = 2 * (map.N2 - num_removed);
                num2 = idx;
                num3 = size(size(A), 2) - map.ii;

                size1 = a_size(1:num1); %ij indices
                size2 = a_size(num1 + 1:num2); %external legs before
                size3 = a_size(num2 + 1:num3); %external legs after
                size4 = a_size(num3 + 1:end); %bond to x0

                perm_vector = [site_ordering_permute(map.N2 - 1); ((2 * map.N2 - 1):size(size(A), 2)).'];
                A = reshape(permute(A, perm_vector), prod(size1), prod(size2), 1, prod(size3), prod(size4));
            end

        case 2

            if map.N - num_removed ~= 0
                perm_vector = [site_ordering_permute(map.N2 - num_removed); ((2 * (map.N2 - num_removed) + 1):size(size(A), 2)).'];

                A = permute(A, perm_vector);

            end

        otherwise
            error("not implemented")
    end

end
