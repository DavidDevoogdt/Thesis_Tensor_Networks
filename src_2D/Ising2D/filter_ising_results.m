function [T_arr, m_arr, marek_arr, corr_arr] = filter_ising_results(T_arr, m_arr, marek_arr, corr_arr, vumps_err_arr, filter, mbound)

    if nargin < 6
        filter = 1;
    end

    if filter == 1
        mask = (T_arr > 0) & (vumps_err_arr < 1e-9) & (vumps_err_arr ~= 0);
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

    if nargin >= 7
        idx2 = find(m_arr < mbound(1));
        idx3 = find(m_arr > mbound(2));
        if isempty(idx2)
            idx2 = numel(idx2);
        end

        if isempty(idx3)
            idx3 = 1;
        end

        T_arr = T_arr(idx3:idx2);
        m_arr = m_arr(idx3:idx2);
        marek_arr = marek_arr(idx3:idx2);
        corr_arr = corr_arr(idx3:idx2);
    end

end
