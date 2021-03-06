function [T_arr, m_arr, marek_arr, corr_arr] = filter_ising_results(T_arr, m_arr, marek_arr, corr_arr, vumps_err_arr, filter)

    if nargin < 6
        filter = 1;
    end

    if filter == 1
        mask = (T_arr > 0) & (vumps_err_arr < 1e-13);
    else
        mask = (T_arr > 0);
    end

    [T_arr, idx] = sort(T_arr(mask));
    m_arr = m_arr(mask);
    m_arr = m_arr(idx);

    marek_arr = marek_arr(mask);
    marek_arr = marek_arr(idx);

    corr_arr = corr_arr(mask);
    corr_arr = corr_arr(idx);
end
