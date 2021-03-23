function [contraction_cell, pat_cells] = get_valid_contractions(obj, map, opts)
    p = inputParser;
    addParameter(p, 'max_index', obj.max_index)
    addParameter(p, 'pattern', {}, @(x) iscell(x)) %additional allowed patterns
    parse(p, opts)

    patterns = p.Results.pattern;
    num_patterns = numel(patterns);

    ie = ndSparse(~cellfun(@isempty, obj.PEPO_cell)) * 1.0;

    iepat = ie;

    for pat_num = 1:num_patterns
        pat = patterns{pat_num};
        ie(pat(1) + 1, pat(2) + 1, pat(3) + 1, pat(4) + 1) = 1;

        iepat(pat(1) + 1, pat(2) + 1, pat(3) + 1, pat(4) + 1) = 0; %remove in case already here
    end

    lookup = cell(size(map.leg_list));
    lookupSize = zeros(map.N, 1);

    con_cell_cell = cell(1, numel(obj.bounds));
    pat_cells_cell = cell(1, numel(obj.bounds));

    total_counter = 0;

    for b = 1:numel(obj.bounds)

        bound = obj.bounds(b);
        tensor_list = cell(1, map.N);
        tensor_list_pat = cell(1, map.N);

        for i = 1:map.N

            T = ie;
            Tpat = iepat;

            connections = map.leg_list{i};

            %only keep sublevel 0 for the given tensors
            if connections(1 + 2) < 0
                T = T(bound, :, :, :);
                Tpat = Tpat(bound, :, :, :);
            end

            if connections(2 + 2) < 0
                T = T(:, bound, :, :);
                Tpat = Tpat(:, bound, :, :);
            end

            if connections(3 + 2) < 0
                T = T(:, :, bound, :);
                Tpat = Tpat(:, :, bound, :);
            end

            if connections(4 + 2) < 0
                T = T(:, :, :, bound);
                Tpat = Tpat(:, :, :, bound);
            end

            A = find(T);
            nA = numel(A);
            lookup{i} = A;
            lookupSize(i) = nA;

            T2 = sparse(1:nA, A, ones(nA, 1), nA, numel(T));
            T2 = ndSparse(T2, [nA, 1, size(T)]);
            tensor_list{i} = T2;

            if nargout > 1
                Apat = find(Tpat);

                c = ismember(A, Apat);

                T2pat = T2;
                if sum(~c) ~= 0
                    T2pat(~c, :) = T2pat(~c, :, :, :, :, :) * 0;
                end

                T2pat = ndSparse(T2pat, [nA, 1, size(T)]);
                tensor_list_pat{i} = T2pat;
            end

        end

        M = ncon_optim(tensor_list, map.leg_list);
        val_con = find(M);

        if nargout > 1
            M_pat = ncon_optim(tensor_list_pat, map.leg_list);
            val_con_pat = find(M_pat);

            pat_cells_cell{b} = ~ismember(val_con, val_con_pat);

        end

        contraction_cell = cell(1, numel(val_con));
        contraction_cell_counter = 1;

        for i = 1:numel(val_con)
            con_cell = cell(1, map.N);

            indices = cell(map.N, 1);
            [indices{:}] = ind2sub(lookupSize, val_con(i));

            for ii = 1:map.N

                real_index = lookup{ii}(indices{ii});

                st = [size(tensor_list{ii}, 3), size(tensor_list{ii}, 4), size(tensor_list{ii}, 5), size(tensor_list{ii}, 6)];

                s = cell(4, 1);

                [s{:}] = ind2sub(st, real_index);

                mask = st == 1;
                aa = find(mask);
                for l = 1:numel(aa)
                    s{aa(l)} = bound;
                end

                con_cell{ii} = [s{:}] - 1; %start indexing at 0
            end

            contraction_cell{contraction_cell_counter} = {con_cell, 'todo'};
            contraction_cell_counter = contraction_cell_counter + 1;
        end

        con_cell_cell{b} = contraction_cell;

        total_counter = total_counter + contraction_cell_counter - 1;
    end

    contraction_cell = cell(total_counter, 1);
    pat_cells = zeros(total_counter, 1) == 1;

    counter = 1;

    for i = 1:numel(obj.bounds)
        num = numel(con_cell_cell{i});
        contraction_cell(counter:counter + num - 1) = con_cell_cell{i};

        m = find(pat_cells_cell{i}) + counter - 1;

        pat_cells(m) = 1;
        counter = counter + num;

    end

    %% old (slow) version

    %     %%%%%%%%% previous results
    %     fixed_mask = p.Results.fixed ~= -1;
    %     num_fixed = sum(fixed_mask);
    %
    %     patterns = p.Results.pattern;
    %     num_patterns = size(patterns, 2);
    %
    %     contraction_cell = cell(0, 1);
    %     contraction_cell_counter = 0;
    %
    %     function [contraction_cell, contraction_cell_counter, vect, correct_index_set] = con_tensors(contraction_cell, contraction_cell_counter, base_level, n, stopind)
    %
    %         tensor_list_indices = cell(1, map.N);
    %
    %         gen_vect = encode_index_array(n, map.internal_legs - num_fixed, stopind) + base_level;
    %         vect = p.Results.fixed;
    %         vect(~fixed_mask) = gen_vect;
    %
    %         correct_index_set = 1;
    %
    %         for i = 1:map.N
    %
    %             if map.is_x_border(i) || map.is_y_border(i)
    %                 legs = [0, 0];
    %
    %                 for j = 1:2
    %                     leg_num = map.leg_list{i}(j);
    %
    %                     if leg_num > 0
    %                         legs(j) = vect(leg_num);
    %                     end
    %                 end
    %
    %                 if map.is_x_border(i)
    %                     O = isempty(obj.boundary_matrix_x{legs(1) + 1, legs(2) + 1});
    %                 else
    %                     O = isempty(obj.boundary_matrix_y{legs(1) + 1, legs(2) + 1});
    %                 end
    %             else
    %
    %                 legs = [0, 0, 0, 0];
    %
    %                 for j = 1:4
    %                     leg_num = map.leg_list{i}(j + 2);
    %
    %                     if leg_num > 0
    %                         legs(j) = vect(leg_num);
    %                     end
    %                 end
    %
    %                 O = isempty(obj.PEPO_cell{legs(1) + 1, legs(2) + 1, legs(3) + 1, legs(4) + 1});
    %
    %                 if O == 1
    %                     for patnum = 1:num_patterns
    %
    %                         if same_pattern(legs, patterns{patnum})
    %                             O = 0;
    %                             break
    %                         end
    %                     end
    %                 end
    %             end
    %
    %             if O == 1
    %
    %                 if obj.visualise == 1
    %                     fprintf("incorrect index set \n");
    %                 end
    %
    %                 correct_index_set = 0;
    %                 break;
    %             end
    %
    %             tensor_list_indices{i} = legs;
    %         end
    %
    %         if correct_index_set
    %
    %             contraction_cell_counter = contraction_cell_counter + 1;
    %             contraction_cell{contraction_cell_counter} = {tensor_list_indices, vect};
    %
    %         end
    %     end
    %
    %     c_index = obj.cycle_index;
    %
    %     tic
    %     if c_index == Inf
    %
    %         for n = 0:(p.Results.max_index + 1)^(map.internal_legs - num_fixed) -1
    %
    %             [contraction_cell, contraction_cell_counter, ~] = con_tensors(contraction_cell, contraction_cell_counter, 0, n, p.Results.max_index);
    %         end
    %
    %     else
    %         fprintf("new");
    %
    %         for n = 0:(obj.cycle_index)^(map.internal_legs - num_fixed) - 1
    %             [contraction_cell, contraction_cell_counter, ~] = con_tensors(contraction_cell, contraction_cell_counter, 0, n, obj.cycle_index - 1);
    %         end
    %
    %         fprintf("cycle_ind");
    %
    %         for n = 0:((p.Results.max_index - obj.cycle_index) +1)^(map.internal_legs - num_fixed) - 1
    %             [contraction_cell, contraction_cell_counter, ~] = con_tensors(contraction_cell, contraction_cell_counter, c_index, n, p.Results.max_index - obj.cycle_index);
    %         end
    %     end
    %     toc

end
