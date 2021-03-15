function [m_arr, T_arr, corr_arr, marek_arr, ctr_arr, vumps_err_arr, pepo_arr] = Ising2D_core(name_prefix, J, g, onsager, chi, i, m_arr, T_arr, corr_arr, marek_arr, ctr_arr, vumps_err_arr, vumps_maxit, pepo_arr)

    if i == 1
         saveboy(sprintf("%s.mat", name_prefix), 'T_arr', 'm_arr', 'corr_arr', 'marek_arr', 'chi', 'J', 'g', 'vumps_err_arr', 'ctr_arr', 'pepo_arr', T_arr, m_arr, corr_arr, marek_arr, chi, J, g, vumps_err_arr, ctr_arr, pepo_arr);
    end

    beta = 1 ./ T_arr(i, :);

    %dir for temp results
    dir_name = sprintf("%s/", name_prefix);
    mkdir(dir_name);

    %hamiltonian setup
    S_x = [0, 1; 1, 0];
    S_y = [0, -1i; 1i, 0];
    S_z = [1, 0; 0, -1];
    I_tensor = eye(2);

    handle = @make_PEPO_2D_A;

    d = 2;
    %J=1;
    %g=0.01;

    H_1_tensor = -J * g * S_x;
    H_2_tensor = -J * (reshape(ncon({S_z, S_z}, {[-1, -3], [-2, -4]}), [d, d, d, d]));

    opts = [];
    opts.testing = 0;
    opts.visualise = 0;
    opts.double = 0;

    % m0 = zeros(size(T0));
    % corr_len0 = zeros(size(T0));
    % marek0 = zeros(size(T0));
    % err0 = zeros(size(T0));
    % ctr0 = zeros(size(T0));

    index = i;

    parfor iii = 1:numel(T_arr(index, :))
        %fprintf("sarted %d:%d:%d T %.4e \n", chi, i, iii, T0(iii));

        T_arr_2 = T_arr(:, iii);
        m_arr_2 = m_arr(:, iii);
        corr_arr_2 = corr_arr(:, iii);
        marek_arr_2 = marek_arr(:, iii);
        ctr_arr_2 = ctr_arr(:, iii);
        vumps_err_arr_2 = vumps_err_arr(:, iii);
        pepo_arr_2 = pepo_arr(:, iii);

        if onsager == 1%for testing purposes

            m_arr_2(index) = m_onsager(T_arr(index, iii), J);
            corr_arr_2(index) = 1;
            marek_arr_2(index) = 1;
            ctr_arr_2(index) = 0;
            vumps_err_arr_2(index) = 1e-13;

            pause(0.1)
        else
            pepo = PEPO(d, -beta(iii) * H_1_tensor, -beta(iii) * H_2_tensor, 5, handle, opts);
            [mm, inv_corr_length, marek_arr_2(index), ctr_arr_2(index), vumps_err_arr_2(index)] = PEPO_get_expectation(pepo, S_z, chi, vumps_maxit, sprintf("%s/vumps_save_v%d_%d_%d.mat", name_prefix, chi, i, iii));

            pepo_arr_2(index) = obj.PEPO_matrix;
            m_arr_2(index) = abs(mm);
            corr_arr_2(index) = 1 / inv_corr_length;

        end

        %T_arr(:, iii) = T_arr_2;
        m_arr(:, iii) = m_arr_2;
        corr_arr(:, iii) = corr_arr_2;
        marek_arr(:, iii) = marek_arr_2;
        ctr_arr(:, iii) = ctr_arr_2;
        vumps_err_arr(:, iii) = vumps_err_arr_2;
        pepo_arr(:, iii) = pepo_arr_2;

        %temp_save in folder
        saveboy(sprintf("%s/temp_%d.mat", name_prefix, iii), 'T_arr_2', 'm_arr_2', 'corr_arr_2', 'marek_arr_2', 'vumps_err_arr_2', 'ctr_arr_2', 'index', 'iii', 'pepo_arr', T_arr_2, m_arr_2, corr_arr_2, marek_arr_2, vumps_err_arr_2, ctr_arr_2, index, iii, pepo_arr_2);
        fprintf("%s %3d:%2d:%2d T:%.4e mag:%.4e xi:%.4e marek gap:%.4f ctr:%3d err:%.4e\n", datestr(now, 'HH:MM:SS'), chi, i, iii, T_arr_2(index), m_arr_2(index), corr_arr_2(index), marek_arr_2(index), ctr_arr_2(index), vumps_err_arr_2(index));

    end

    %save to main variables and delete temp dir
    saveboy(sprintf("%s.mat", name_prefix), 'T_arr', 'm_arr', 'corr_arr', 'marek_arr', 'chi', 'J', 'g', 'vumps_err_arr', 'ctr_arr', 'pepo_arr', T_arr, m_arr, corr_arr, marek_arr, chi, J, g, vumps_err_arr, ctr_arr, pepo_arr);

    delete(sprintf("%s/temp_*.mat", name_prefix));
    %rmdir(dir_name);

end

function m = m_onsager(T, J)
    T_c = 2 * J / (log(1 + sqrt(2)));
    m = T;
    mask = T < T_c;
    m(mask) = (1 - sinh((2 * J) ./ T(mask)).^(-4)).^(1/8);
    m(~mask) = 0;
end
