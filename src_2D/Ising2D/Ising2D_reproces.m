function Ising2D_reproces(suffix)

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');
    pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
    fold2 = strjoin(pathparts, '/');

    filenames = {
            'Ising2D_g=2.5000e+00_chi=8_02_April_2021_09:38';
            'Ising2D_g=2.5000e+00_chi=11_02_April_2021_09:38';
            'Ising2D_g=2.5000e+00_chi=16_02_April_2021_09:38';
            'Ising2D_g=2.5000e+00_chi=23_02_April_2021_09:38';
            'Ising2D_g=2.5000e+00_chi=32_02_April_2021_09:38';
            'Ising2D_g=2.5000e+00_chi=45_02_April_2021_09:38';
            'Ising2D_g=1.5000e+00_chi=8_02_April_2021_11:17';
            'Ising2D_g=1.5000e+00_chi=11_02_April_2021_11:17';
            'Ising2D_g=1.5000e+00_chi=16_02_April_2021_11:17';
            'Ising2D_g=1.5000e+00_chi=23_02_April_2021_11:17';
            'Ising2D_g=1.5000e+00_chi=32_02_April_2021_11:17';
            'Ising2D_g=1.5000e+00_chi=45_02_April_2021_11:17';
            'Ising2D_g=1.5000e+00_chi=64_02_April_2021_11:17';
            'Ising2D_g=0.0000e+00_chi=8_02_April_2021_11:46';
            'Ising2D_g=0.0000e+00_chi=11_02_April_2021_11:46';
            'Ising2D_g=0.0000e+00_chi=16_02_April_2021_11:46';
            'Ising2D_g=0.0000e+00_chi=23_02_April_2021_11:46';
            'Ising2D_g=0.0000e+00_chi=32_02_April_2021_11:46';
            'Ising2D_g=0.0000e+00_chi=45_02_April_2021_11:46';
    %'Ising2D_g=2.5000e+00_chi=20_22_March_2021_10:11';
            };

    reproces_opts = struct('doEpsi', 0, 'doVumps', 0);

    dt = datestr(now, 'dd_mmmm_yyyy_HH:MM');

    parfor i = 1:numel(filenames)

        %load template
        [~, template] = fetch_matfiles(filenames{i}, struct);

        if template.found == 0 %guess the settings, legacy
            vumps_opts = [];
            vumps_opts.vumps_maxit = 1300;
            vumps_opts.tolfixed = 1e-10;

            template.vumps_opts = vumps_opts;

            template.handle = @make_PEPO_2D_A;
            S_z = [1, 0; 0, -1];
            template.X = S_z; %observable

            template.model_params = models('t_ising', struct('g', 2.5));
            template.pepo_opts = struct();
        end

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

        fetch_matfiles(filenames{i}, opts)

    end
end

function call_back_fn(results, save_vars, baseFileName, template, reproces_opts)

    if ~isfield(template.vumps_opts, 'chi_max')
        template.vumps_opts.chi_max = save_vars.A{1}.dims(1);
    end

    save_vars.fname = strrep (strrep(baseFileName, 'save_vars_', ''), '.mat', '');

    good_point = (results.T > 0) & (results.err < template.vumps_opts.tolfixed) & (results.err ~= 0) & (results.inv_corr_length > 1e-5);

    if good_point
        try
            Ising2D_core(save_vars, template, results.T, results, reproces_opts)
        catch
            %fprintf("faild: %s %s \n", template.name_prefix, baseFileName)
        end
    end
end
