function sampling_reproces(suffix)
    %function to redo calculations

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');
    pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
    fold2 = strjoin(pathparts, '/');

    filenames = {
            'TIM_T=0.7_order_5_chi=8_trunc_20_sym=1_18_May_2021_19:41';
            };

    reproces_opts = struct('doEpsi', 0, 'doVumps', 0);
    dt = datestr(now, 'dd_mmmm_yyyy_HH:MM');

    for i = 1:numel(filenames)

        %load template
        [~, template] = sampling_fetch(filenames{i}, struct);

        if strcmpi(suffix, 'none') %update inplace
            template.name_prefix = sprintf("%s/%s", fold2, filenames{i});
        else
            template.name_prefix = sprintf("%s/%s_%s", fold2, filenames{i}, suffix);
        end

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

        sampling_fetch(filenames{i}, opts)

    end
end

function call_back_fn(results, save_vars, baseFileName, template, reproces_opts)
    save_vars.fname = strrep (strrep(baseFileName, 'save_vars_', ''), '.mat', '');
    good_point = (results.T > 0) & (results.err < template.vumps_opts.tolfixed) & (results.err ~= 0) & (results.inv_corr_length > 1e-5);
    if good_point
        sampling_core(save_vars, template, results.T, results, reproces_opts)
    end
end
