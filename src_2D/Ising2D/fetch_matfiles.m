function [T_arr, m_arr, marek_arr, corr_arr, vumps_err_arr, ctr_arr, J, chi, g,pepo_arr] = fetch_matfiles(name_prefix)

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');
    pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
    fold2 = strjoin(pathparts, '/');

    name = sprintf("%s.mat", name_prefix);
    folder = sprintf("%s/%s/", fold2, name_prefix);
    
    if nargout <10
        load(name, 'T_arr', 'm_arr', 'corr_arr', 'marek_arr', 'chi', 'J', 'g', 'vumps_err_arr', 'ctr_arr');
    else
        load(name, 'T_arr', 'm_arr', 'corr_arr', 'marek_arr', 'chi', 'J', 'g', 'vumps_err_arr', 'ctr_arr','pepo_arr');
    end

    if exist(folder, 'dir')
        myFiles = dir(fullfile(folder, 'temp_*.mat')); %gets all wav files in struct
        for k = 1:length(myFiles)
            baseFileName = myFiles(k).name;
            fullFileName = fullfile(folder, baseFileName);

            load(fullFileName, 'T_arr_2', 'm_arr_2', 'corr_arr_2', 'marek_arr_2', 'vumps_err_arr_2', 'ctr_arr_2', 'index', 'iii');

            mask = T_arr(:, iii) == 0;

            T_arr(:, iii) = T_arr_2;
            m_arr(:, iii) = m_arr_2;
            corr_arr(:, iii) = corr_arr_2;
            marek_arr(:, iii) = marek_arr_2;
            ctr_arr(:, iii) = ctr_arr_2;
            vumps_err_arr(:, iii) = vumps_err_arr_2;
        end
    end
% 
%     T_arr = reshape(T_arr, [], 1);
%     m_arr = reshape(m_arr, [], 1);
%     marek_arr = reshape(marek_arr, [], 1);
%     corr_arr = reshape(corr_arr, [], 1);
%     vumps_err_arr = reshape(vumps_err_arr, [], 1);
%     ctr_arr = reshape(ctr_arr, [], 1);

end
