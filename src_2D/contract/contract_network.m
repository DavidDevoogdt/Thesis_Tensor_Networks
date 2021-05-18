function M = contract_network(obj, map, opts)
    %contracts network for a given geometry. matrix determines wether is is contracted at once or cell per cell
    %results can be traced before contraction.

    p = inputParser;
    addParameter(p, 'matrix', 0)
    addParameter(p, 'lnprefact', obj.nf);
    addParameter(p, 'trace', false);
    addParameter(p, 'non_trace_num', 0);
    parse(p, opts)

    if p.Results.matrix == 0

        M = zeros(dimension_vector(obj.dim, 2 * map.N2));
        con_cells = get_valid_contractions(obj, map, struct);
        M = -contract_con_cells(obj, map, p.Results.lnprefact, M, con_cells, struct('permute', 0));

    else
        tensor_list = cell(1, map.N);

        mult_fact = exp(p.Results.lnprefact - obj.nf);

        T0 = obj.PEPO_matrix / mult_fact;

        for i = 1:map.N
            T = T0;
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

            if p.Results.trace == true && i ~= p.Results.non_trace_num
                T = ncon({T}, {[1, 1, -3, -4, -5, -6, -1, -2]}); %trace out physical indices
            end

            tensor_list{i} = T;
        end

        if map.N < 6
            M = ncon(tensor_list, map.leg_list);
        else
            M = ncon_optim(tensor_list, map.leg_list);
        end

    end

end
