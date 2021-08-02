function data = sampling_filter(data, opts)
    %filter out bad points and return the data

    if ~isfield(opts, 'tol')
        opts.tol = 1e-10;
    end

    if ~isfield(opts, 'Tbound')
        opts.Tbound = [0, Inf];
    end

    free_var = data.free_var;

    mask = (data.(free_var) > opts.Tbound(1)) & (data.(free_var) < opts.Tbound(2)) & (data.err < opts.tol) & (data.err ~= 0);

    [~, idx] = sort(data.(free_var)(mask));

    for i = 1:numel(data.fields)
        nf = data.fields{i};
        data.(nf) = data.(nf)(mask, :);
        data.(nf) = data.(nf)(idx, :);
    end
end
