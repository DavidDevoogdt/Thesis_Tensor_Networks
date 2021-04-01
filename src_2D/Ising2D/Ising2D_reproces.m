function Ising2D_reproces

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

    for i = 1:numel(filenames)

        root_folder = sprintf("%s/%s", fold2, filenames{i});
        fprintf("root fold: %s \n", root_folder);

        data = fetch_matfiles(filenames{i}, struct);

        for j = 1:numel(s)

            fprintf("%s \n", this_file);

            if ~exist(this_file, 'file')
                T_arr(m, n) = 0;
                warning('%s missing', this_file);
            else

                S = load(this_file, 'G0', 'A');
                %
                S_z = [1, 0; 0, -1];
                %
                %                 try
                %                     [mm, inv_corr_length, delta] = PEPO_get_expectation ([], S_z, [], [], [], S.A, S.G0, pepo_arr{m, n});
                %                 catch
                %                     T_arr(m, n) = 0;
                %                     warning('%s calculating environment failed,removing dat point', this_file);
                %                     continue;
                %                 end
                %
                %                 %m_arr(s(j)) = abs(mm);
                %                 %corr_arr(s(j)) = 1 / inv_corr_length;
                %                 %m_arr(s(j) = delta;

                template.J = J;
                template.g = g;
                template.chi = chi;
                template.vumps_maxit = -1;
                template.X = S_z;

                save_vars.A = S.A;
                save_vars.G0 = S.G0;
                save_vars.PEPO_matrix = pepo_arr{m, n};

                results.T = T_arr(s(j));
                results.m = m_arr(s(j));
                results.marek = marek_arr(s(j));
                %results.eps_i = [];
                results.ctr = ctr_arr(s(j));
                results.err = vumps_err_arr(s(j));
                results.inv_corr_length = 1 ./ corr_arr(s(j));

                saveboy(sprintf("%s/results_%d_%d.mat", root_folder, m, n), 'results', 'template', results, template);
                saveboy(sprintf("%s/save_vars_%d_%d.mat", root_folder, m, n), 'save_vars', save_vars);

            end

        end
        %         %save again
        %         name_rep = sprintf("%s/%s_rep.mat", fold2, filenames{i});
        %
        %         saveboy(name_rep, 'T_arr', 'm_arr', 'corr_arr', 'marek_arr', 'chi', 'J', 'g', 'vumps_err_arr', 'ctr_arr', T_arr, m_arr, corr_arr, marek_arr, chi, J, g, vumps_err_arr, ctr_arr);
        %         fprintf('done with %s\n', root_folder);
    end
end
