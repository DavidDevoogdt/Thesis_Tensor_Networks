function M = contract_network(obj, map, opts)
    %generate all index sets for a given configuration

    %             struct('max_index', [] ,'matrix', [] ,'fixed', []);

    p = inputParser;
    addParameter(p, 'max_index', obj.max_index)
    addParameter(p, 'matrix', 0)
    addParameter(p, 'fixed', zeros(map.internal_legs, 1) - 1)
    addParameter(p, 'lnprefact', obj.nf);
    parse(p, opts)

    M = zeros(dimension_vector(obj.dim, 2 * map.N2));
    %M2 = zeros(dimension_vector(obj.dim, 2 * map.N2));

    if p.Results.matrix == 0
        %         correct_index_sets = get_valid_contractions(obj, map, opts);
        %
        %
        %         for i = 1:numel(correct_index_sets)
        %             iset = correct_index_sets{i};
        %             %vect = iset{2};
        %             legs = iset{1};
        %
        %             tensors = fetch_PEPO_cells(obj, map, legs,p.Results.lnprefact);
        %
        %             M = M + ncon(tensors, map.leg_list);
        %         end

        con_cells = get_valid_contractions(obj, map, opts);
        M = -contract_con_cells(obj, map, p.Results.lnprefact, M, con_cells, struct('permute', 0));

    else
        tensor_list = cell(1, map.N);

        mult_fact = exp(p.Results.lnprefact - obj.nf);

        for i = 1:map.N
            T = obj.PEPO_matrix / mult_fact;
            connections = map.leg_list{i};
            %only keep sublevel 0 for the given tensors
            if connections(1 + 2) < 0
                T = T(:, :, 1, :, :, :);
            end

            if connections(2 + 2) < 0
                T = T(:, :, :, 1, :, :);
            end

            if connections(3 + 2) < 0
                T = T(:, :, :, :, 1, :);
            end

            if connections(4 + 2) < 0
                T = T(:, :, :, :, :, 1);
            end

            tensor_list{i} = T;
        end

        M = ncon_optim(tensor_list, map.leg_list);

    end

end
