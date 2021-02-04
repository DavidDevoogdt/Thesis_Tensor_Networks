function contraction_cell = get_valid_contractions(obj, map, opts)
    p = inputParser;
    addParameter(p, 'max_index', obj.max_index)
    %addParameter(p, 'matrix', 0)
    %addParameter(p, 'fixed', zeros(map.internal_legs, 1) - 1)
    addParameter(p, 'pattern', {}, @(x) iscell(x))%additional allowed patterns
    parse(p, opts)

    tensor_list = cell(1, map.N);
    patterns = p.Results.pattern;
    num_patterns = size(patterns, 2);

    ie = ndSparse(~cellfun(@isempty, obj.PEPO_cell));

    for pat_num = 1:num_patterns
        pat = patterns{pat_num};
        ie(pat(1) + 1, pat(2) + 1, pat(3) + 1, pat(4) + 1) = 1;
    end

    lookup = cell(size(map.leg_list));
    lookupSize = zeros(map.N, 1);

    for i = 1:map.N

        T = ie;
        connections = map.leg_list{i};

        %only keep sublevel 0 for the given tensors
        if connections(1 + 2) < 0
            T = T(1, :, :, :);
        end

        if connections(2 + 2) < 0
            T = T(:, 1, :, :);
        end

        if connections(3 + 2) < 0
            T = T(:, :, 1, :);
        end

        if connections(4 + 2) < 0
            T = T(:, :, :, 1);
        end

        A = find(T);
        nA = numel(A);
        lookup{i} = A;
        lookupSize(i) = nA;

        T2 = sparse(1:nA, A, ones(nA, 1), nA, numel(T));
        T2 = ndSparse(T2, [nA, 1, size(T)]);
        tensor_list{i} = T2;
    end

    M = ncon(tensor_list, map.leg_list);
    val_con = find(M);

    contraction_cell = cell(1, numel(val_con));
    contraction_cell_counter = 1;

    for i = 1:numel(val_con)
        con_cell = cell(1, map.N);

        indices = cell(map.N, 1);
        [indices{:}] = ind2sub(lookupSize, val_con(i));

        for ii = 1:map.N

            real_index = lookup{ii}(indices{ii});

            st = size(tensor_list{ii});

            [s1, s2, s3, s4] = ind2sub(st(3:end), real_index);
            con_cell{ii} = [s1, s2, s3, s4] - 1; %start indexing at 0
        end

        contraction_cell{contraction_cell_counter} = {con_cell, 'todo'};
        contraction_cell_counter = contraction_cell_counter + 1;
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
