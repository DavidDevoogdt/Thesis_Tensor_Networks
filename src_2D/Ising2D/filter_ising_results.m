function data = filter_ising_results(data, opts)

    if ~isfield(opts, 'tol')
        opts.tol = 1e-10;
    end

    mask = (data.T > 0) & (data.err < opts.tol) & (data.err ~= 0);

    [~, idx] = sort(data.T(mask));

    for i = 1:numel(data.fields)
        nf = data.fields{i};
        data.(nf) = data.(nf)(mask, :);
        data.(nf) = data.(nf)(idx, :);
    end
end
