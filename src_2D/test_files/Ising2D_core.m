function [m0, corr_len0, marek0, T_max] = Ising2D_core(T0, J, g, onsager, chi, i, T_max)

    beta = 1 ./ T0;

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

    m0 = zeros(size(T0));
    corr_len0 = zeros(size(T0));
    marek0 = zeros(size(T0));
     

    T_max = ones(size(T0)) .* T_max;

    parfor iii = 1:numel(T0)
        %fprintf("sarted %d:%d:%d T %.4e \n", chi, i, iii, T0(iii));

        ctr = 0;
        err = 0;
        if onsager == 1
            m0(iii) = m_onsager(T0(iii), J);
            corr_len0(iii) = 1;
            marek0(iii) = 1;

            if m0(iii)<0.2
                err = 1;
            end
            
            pause(0.1)
        else

            pepo = PEPO(d, -beta(iii) * H_1_tensor, -beta(iii) * H_2_tensor, 5, handle, opts);
            [mm, inv_corr_length, marek0(iii), ctr, err] = PEPO_get_expectation(pepo, S_z, chi);
            m0(iii) = abs(mm);
            corr_len0(iii) = 1 / inv_corr_length;

           
        end
        
%          if err == 1
%                 T_max(iii) = 1 / beta(iii);
%                 m0(iii) = 1e-15;
%                 corr_len0(iii) = -1;
%                 marek0(iii)=-1;
%          end

        
        fprintf("%d:%d:%d T %.4e mag:%.4e xi:%.4e marek gap:%.4f, ctr %d,err %d \n", chi, i, iii, T0(iii), m0(iii), corr_len0(iii), marek0(iii), ctr, err);

    end

    T_max = min(T_max);

end

function m = m_onsager(T, J)
    T_c = 2 * J / (log(1 + sqrt(2)));
    m = T;
    mask = T < T_c;
    m(mask) = (1 - sinh((2 * J) ./ T(mask)).^(-4)).^(1/8);
    m(~mask) = 0;
end
