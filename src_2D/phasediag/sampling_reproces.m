function sampling_reproces(name, suffix)
    %function to redo calculations

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');
    pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
    fold2 = strjoin(pathparts, '/');

    if ~ischar(name)
        nn = Ising2D_names(name);
        filenames = nn{1};
    else
        filenames = {name};
    end

    reproces_opts = struct; %('doEpsi', 1, 'doVumps', 0);
    dt = datestr(now, 'dd_mmmm_yyyy_HH:MM');

    for i = 1:numel(filenames)

        %load template
        [~, template] = sampling_fetch(filenames{i}, struct);

        if strcmpi(suffix, 'none') %update inplace
            template.name_prefix = sprintf("%s/%s", fold2, filenames{i});
        else
            template.name_prefix = sprintf("%s/%s_%s", fold2, filenames{i}, suffix);
        end

        template.vumps_opts.tolfixed = 1e-9;
        template.vumps_opts.vumps_maxit = 5;
        
        
        dir_name = sprintf("%s/", template.name_prefix);
        if ~exist(dir_name, 'dir')
            mkdir(dir_name);
        end

        saveboy(sprintf("%s/template.mat", template.name_prefix), 'template', template);
        fprintf("starting with %s \n", template.name_prefix);
        template.time = dt;

        opts = [];
        opts.call_back_fn = @(y, z, a) call_back_fn(y, z, a, template, reproces_opts);
        opts.save_vars = 1;

        sampling_fetch(filenames{i}, opts);

    end
end

function call_back_fn(results, save_vars, baseFileName, template, reproces_opts)
    save_vars.fname = strrep (strrep(baseFileName, 'save_vars_', ''), '.mat', '');
    good_point = (results.T > 0) & (results.err < 1e-10  ) & (results.err ~= 0);
    
    if good_point
        
        v = sampling_core(save_vars, template, results.T, results, reproces_opts);
        fprintf("%s",v);
    end
end
