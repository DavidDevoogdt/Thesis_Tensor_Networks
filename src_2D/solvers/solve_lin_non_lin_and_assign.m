function obj = solve_lin_non_lin_and_assign(obj, map, patterns, ln_prefact, opts)

    obj = fill_rand(obj, patterns);

    p = inputParser;
    addParameter(p, 'maxit', 20)
    addParameter(p, 'solved', 1e-13)
    addParameter(p, 'minchange', 0.99)
    addParameter(p, 'display', 0)
    parse(p, opts);

    maxit = p.Results.maxit;
    
    prev_err =  calculate_error(obj, map.num_map, obj.numopts);

    for j = 1:maxit

        for i = 1:numel(patterns)
            obj = solve_lin_and_assign(obj, map, patterns(i), ln_prefact,-1,1);
        end
        
        err = calculate_error(obj, map.num_map, obj.numopts);
        
        if p.Results.display==1
           fprintf("%3d: %.4e\n",j,err); 
        end
        
        if err<p.Results.solved
           break; 
        end
        
        if err/prev_err > p.Results.minchange
            break;
        end
        
        prev_err = err;
    end
end
