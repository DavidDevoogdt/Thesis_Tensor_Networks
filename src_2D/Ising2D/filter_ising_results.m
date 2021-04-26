function data = filter_ising_results(data, opts)

    if ~isfield(opts, 'tol')
        opts.tol = 1e-10;
    end

    if ~isfield(opts, 'Tbound')
        opts.Tbound = [0,Inf];
    end

    mask = (data.T > opts.Tbound(1)) & ( data.T < opts.Tbound(2) ) & (data.err < opts.tol) & (data.err ~= 0)  ;
    


    [~, idx] = sort(data. (data.free_var)(mask));

    for i = 1:numel(data.fields)
        nf = data.fields{i};
        data.(nf) = data.(nf)(mask, :);
        data.(nf) = data.(nf)(idx, :);
    end
end
