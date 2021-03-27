function Ising2D_reprocess

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');
    pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
    fold2 = strjoin(pathparts, '/');

    filenames = {
            'Ising2D_g=2.5000e+00_chi=20_22_March_2021_10:11';
            'Ising2D_g=2.5000e+00_chi=25_22_March_2021_10:11';
            'Ising2D_g=2.5000e+00_chi=30_22_March_2021_10:11';
            'Ising2D_g=2.5000e+00_chi=40_22_March_2021_10:46';
            'Ising2D_g=2.5000e+00_chi=60_22_March_2021_11:49';
            'Ising2D_g=2.5000e+00_chi=65_22_March_2021_13:52';
            'Ising2D_g=2.5000e+00_chi=70_22_March_2021_13:51';
            };

    parfor i = 1:numel(filenames)

        root_folder = sprintf("%s/%s", fold2, filenames{i});
        fprintf("root fold: %s \n", root_folder);

        [T_arr, m_arr_orig, ~, ~, vumps_err_arr, ctr_arr, J, chi, g, pepo_arr] = fetch_matfiles(filenames{i});

        m_arr = zeros(size(T_arr));
        marek_arr = zeros(size(T_arr));
        corr_arr = zeros(size(T_arr));

        if ~exist(root_folder, 'dir')
            error('vumps_save not found');
        end

        s = find(m_arr_orig ~= 0);

        for j = 1:numel(s)

            [m, n] = ind2sub(size(T_arr), s(j));

            this_file = sprintf("%s/vumps_save_v%d_%d_%d.mat", root_folder, chi, m, n);

            fprintf("%s \n", this_file);

            if ~exist(this_file, 'file')
                T_arr(m, n) = 0;
                warning('missing file');
            else

                S = load(this_file, 'G0', 'A');

                S_z = [1, 0; 0, -1];

                [mm, inv_corr_length, delta] = PEPO_get_expectation ([], S_z, [], [], [], S.A, S.G0, pepo_arr{m, n});

                m_arr(s(j)) = abs(mm);
                corr_arr(s(j)) = 1 / inv_corr_length;
                marek_arr(s(j)) = delta;
            end
        end
        %save again
        name_rep = sprintf("%s/%s_rep.mat", fold2, filenames{i});

        saveboy(name_rep, 'T_arr', 'm_arr', 'corr_arr', 'marek_arr', 'chi', 'J', 'g', 'vumps_err_arr', 'ctr_arr', T_arr, m_arr, corr_arr, marek_arr, chi, J, g, vumps_err_arr, ctr_arr);
        fprintf('done with %s\n', root_folder);
    end
end
