function target2 = contract_con_cells(obj, map, ln_prefactor, target2, con_cells, workers)
    if nargin < 6
        workers = 8;
    end

    n = numel(con_cells);

    sizes = zeros(workers, 1);
    sizes(1:workers - 1) = floor(n / workers);
    sizes(workers) = n - sum(sizes);
    start_vals = [0; cumsum(sizes(1:end - 1))];

    targets = cell(workers, 1);

    parfor j = 1:workers
        target = zeros(size(target2));

        start_val = start_vals(j);

        for i = 1:sizes(j)
            legs = con_cells{i + start_val}{1};

            temp_list_1 = fetch_PEPO_cells(obj, map, legs, ln_prefactor);

            A1 = ncon(temp_list_1, map.leg_list);
            perm_vect = [site_ordering_permute(map.N2); (2 * map.N2 + 1:size(size(A1), 2))'];

            target = target - reshape(permute(A1, perm_vect), size(target));

        end

        targets{j} = target;
    end

    for j = 1:workers
        target2 = target2 + targets{j};
    end

end
