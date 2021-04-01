function Ising2D_reproces(suffix)

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');
    pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
    fold2 = strjoin(pathparts, '/');

    filenames = {
            'Ising2D_g=2.5000e+00_chi=20_22_March_2021_10:11';
    % 'Ising2D_g=2.5000e+00_chi=25_22_March_2021_10:11';
    % 'Ising2D_g=2.5000e+00_chi=30_22_March_2021_10:11';
    % 'Ising2D_g=2.5000e+00_chi=40_22_March_2021_10:46';
    % 'Ising2D_g=2.5000e+00_chi=60_22_March_2021_11:49';
    % 'Ising2D_g=2.5000e+00_chi=65_22_March_2021_13:52';
    % 'Ising2D_g=2.5000e+00_chi=70_22_March_2021_13:51';
            };

    dt = datestr(now, 'dd_mmmm_yyyy_HH:MM');

    for i = 1:numel(filenames)

        settings.name_prefix = sprintf("%s/%s_%s", fold2, filenames{i}, suffix);
        settings.time = dt;

        opts.call_back_fn = @(x, y, z, a) call_back_fn(x, y, z, a, settings);
        opts.save_vars = 1;

        fetch_matfiles(filenames{i}, opts)

    end
end

function call_back_fn(template, results, save_vars, baseFileName, settings)

    template.fname = baseFileName;

    vumps_opts.chi_max = save_vars.A{1}.dims(1);
    vumps_opts.vumps_maxit = 1300;
    vumps_opts.tolfixed = 1e-10;

    template.vumps_opts = vumps_opts;
    template.handle = @make_PEPO_2D_A;

    template.name_prefix = settings.name_prefix;

    save_vars.fname = strrep (strrep(baseFileName, 'save_vars_', ''), '.mat', '');

    template.time = settings.time;

    S_z = [1, 0; 0, -1];
    template.X = S_z; %observable

    template.model_params = models('t_ising', struct('g', 2.5));
    template.pepo_opts = struct();

    Ising2D_core(save_vars, template, results.T, results)

end
