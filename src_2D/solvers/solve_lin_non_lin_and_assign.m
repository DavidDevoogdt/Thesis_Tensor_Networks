function obj = solve_lin_non_lin_and_assign(obj, map, patterns, ln_prefact, opts, patfilt )

    

    obj = fill_rand(obj, patterns);

    p = inputParser;
    addParameter(p, 'maxit', 20)
    addParameter(p, 'solved', 1e-13)
    addParameter(p, 'minchange', 0.99)
    addParameter(p, 'display', 0)
    parse(p, opts);

    maxit = p.Results.maxit;

    prev_err = calculate_error(obj, map.num_map, obj.numopts);

    [all_con_cells, pat_cells] = get_valid_contractions(obj, map, struct('max_index', obj.current_max_index, 'pattern', {patterns}));
    
    
    
    if nargin>=6
        pcells = find(pat_cells);
        
        matching = zeros(size(pcells));
        for i=1:numel(pcells)
            ccell = all_con_cells{pcells(i)}{1};
            
            match = zeros(size(ccell));
            for j=1:numel(ccell)
                for k=1:numel(patfilt)
                    if sum(patfilt{k}==ccell{j}) == 4
                        match(j)=1;
                        break;
                    end
                end
            end
            if match == ones(size(match))
               matching(i)=1; 
            end
        end  
        
        pat_cells( pcells( ~matching ) ) = 0;
        
    end
    
    for j = 1:maxit

        for i = 1:numel(patterns)
            obj = solve_lin_and_assign(obj, map, patterns(i), ln_prefact, -1, 1,all_con_cells, pat_cells);
        end

        err = calculate_error(obj, map.num_map, obj.numopts);

        if p.Results.display == 1
            fprintf("%3d: %.4e\n", j, err);
        end

        if err < p.Results.solved
            break;
        end

        if err / prev_err > p.Results.minchange
            break;
        end

        prev_err = err;
    end  
       

end
