function [contraction_cell, pat_cells] = get_valid_contractions(obj, map, opts)
    %calculates all possible combinations of the virtual levels. Under the hood, the problem is reprocessed into a PEPS network, contracted and read out.
    % addition patterns are added to possible contractions. if 2 output arguments are requested, an array indexing all contractions involving these patterns is given

    p = inputParser;
    addParameter(p, 'pattern', {}, @(x) iscell(x)) %additional allowed patterns
    parse(p, opts)

    patterns = p.Results.pattern;
    num_patterns = numel(patterns);

    ie = ndSparse(~cellfun(@isempty, obj.PEPO_cell)) * 1.0;

    iepat = ie;

    for pat_num = 1:num_patterns
        pat = patterns{pat_num};
        ie(pat(1) + 1, pat(2) + 1, pat(3) + 1, pat(4) + 1) = 1;
        iepat(pat(1) + 1, pat(2) + 1, pat(3) + 1, pat(4) + 1) = 0;
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

        %make PEPS network
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

        %contract PEPS network
        if map.N > 6
            M = ncon_optim(tensor_list, map.leg_list);
        else
            M = ncon(tensor_list, map.leg_list);
        end

        val_con = find(M);

        %contract PEPS network without extra patterns and find difference
        if nargout > 1
            if map.N > 6
                M_pat = ncon_optim(tensor_list_pat, map.leg_list);
            else
                M_pat = ncon(tensor_list_pat, map.leg_list);
            end

            val_con_pat = find(M_pat);

            pat_cells_cell{b} = ~ismember(val_con, val_con_pat);

        end

        %precompute size to gain speed
        st2 = ones(map.N, 4);
        for ii = 1:map.N
            sz = size(tensor_list{ii});
            nsz = numel(sz) - 2;

            st2(ii, 1:nsz) = sz(3:end);
        end

        contraction_cell = cell(1, numel(val_con));
        contraction_cell_counter = 1;

        %decode the values to the original indices
        for i = 1:numel(val_con)
            con_cell = cell(1, map.N);

            indices = cell(map.N, 1);
            [indices{:}] = ind2sub(lookupSize, val_con(i));

            for ii = 1:map.N
                real_index = lookup{ii}(indices{ii});

                st = st2(ii, :);
                s = cell(4, 1);

                [s{:}] = ind2sub(st, real_index);

                mask = st == 1;
                aa = find(mask);
                for l = 1:numel(aa)
                    s{aa(l)} = bound;
                end

                con_cell{ii} = [s{:}] - 1;
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
end
