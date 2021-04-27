function target2 = contract_con_cells(obj, map, ln_prefactor, target2, con_cells, opts)
    if nargin < 6
        opts = struct();
    end

    p = inputParser;
    addParameter(p, 'workers', 1)
    addParameter(p, 'permute', 1)
    parse(p, opts)

    workers = p.Results.workers;
    perm = p.Results.permute;

    n = numel(con_cells);

    sizes = zeros(workers, 1);
    sizes(1:workers - 1) = floor(n / workers);
    sizes(workers) = n - sum(sizes);
    start_vals = [0; cumsum(sizes(1:end - 1))];

    targets = cell(workers, 1);

    cc = cell(workers, 1);

    for j = 1:workers
        start_val = start_vals(j);
        cc{j} = con_cells(start_val + 1:start_val + sizes(j));
    end

    for j = 1:workers
        target = zeros(size(target2));

        %start_val = start_vals(j);

        c = cc{j};

        for i = 1:sizes(j)
            legs = c{i}{1};

            temp_list_1 = fetch_PEPO_cells(obj, map, legs, ln_prefactor);

            if sum(cellfun(@isempty, temp_list_1)) ==0 %uninitialised cell, skip
            
                A1 = ncon(temp_list_1, map.leg_list);

                if perm == 1
                    perm_vect = [site_ordering_permute(map.N2); (2 * map.N2 + 1:size(size(A1), 2))'];
                    target = target - reshape(permute(A1, perm_vect), size(target));
                else
                    target = target - A1;
                end
            end

        end

        targets{j} = target;
    end

    for j = 1:workers
        target2 = target2 + targets{j};
    end

end
