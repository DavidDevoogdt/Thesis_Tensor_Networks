%for tensors:       ( beta)
%            (alpha)-- O -- (gamma)  = O(i,j,alpha,beta,gamma,delta)
%                    (delta)             1 2   3     4    5     6
% for tensors containig multiple ij's: O numbered from left to right and
% for a given vertical pos from up till down

classdef PEPO

    properties
        dim
        H_1_tensor
        H_2_tensor
        PEPO_cell
        nf %normalisation factor
        max_index
        testing
        visualise
        virtual_level_sizes_horiz
        virtual_level_sizes_vert
        PEPO_matrix
        current_max_index
        numopts
        cycleopts
        cycle_index
        boundary_matrix_x
        boundary_matrix_y
        boundary_vect
        bounds
        order
        inv_eps
        complex
        err_tol
        error_code
    end

    methods

        function [obj, err_code] = PEPO(d, H_1_tensor, H_2_tensor, order, make_PEPO_handle, opts)
            numopts.numbered = 1;
            obj.numopts = numopts;

            cycleopts.numbered = 1;
            cycleopts.v_cyclic = 0;
            cycleopts.h_cyclic = 1;
            obj.cycleopts = cycleopts;

            obj.dim = d;
            obj.H_1_tensor = H_1_tensor;
            obj.H_2_tensor = H_2_tensor;

            obj.complex = false;

            %parse opts
            p = inputParser;
            addParameter(p, 'testing', 0)
            addParameter(p, 'visualise', 0)
            addParameter(p, 'double', 0)
            addParameter(p, 'inv_eps', 1e-12)
            addParameter(p, 'err_tol', 1e-13)

            parse(p, opts)

            obj.err_tol = p.Results.err_tol;

            obj.inv_eps = p.Results.inv_eps;

            if mod(order, 2)
                max_index = (order - 1) / 2;
            else
                max_index = order / 2;
            end

            if p.Results.double == 1
                max_index = 2 * max_index + 1;
            end

            obj.PEPO_cell = cell(max_index + 8, max_index + 8, max_index + 8, max_index + 8);
            obj.boundary_matrix_x = cell(max_index + 1, max_index + 2);
            obj.boundary_matrix_y = cell(max_index + 1, max_index + 2);
            obj.boundary_matrix_x{1, 1} = reshape(1, 1, 1);
            obj.boundary_matrix_y{1, 1} = reshape(1, 1, 1);

            obj.virtual_level_sizes_horiz = 1;
            obj.virtual_level_sizes_vert = 1;

            obj.max_index = max_index;

            obj.order = order;
            obj.cycle_index = Inf;

            obj.testing = p.Results.testing;
            obj.visualise = p.Results.visualise;

            %calculate ln of normailisation fac
            map = create_map(1:2, numopts);
            [~, nf2] = H_exp(obj, map, 0, true);
            obj.nf = nf2;

            obj.PEPO_cell{1, 1, 1, 1} = reshape(eye(d) / exp(obj.nf), [d, d, 1, 1, 1, 1]);

            %non generic PEPO code

            [obj, err_code] = make_PEPO_handle(obj);

            obj = cell2matrix(obj); %save matrix form

            obj.error_code = err_code;

            if err_code == 1
                warning("PEPO costruction failed, try lower beta");
            end

        end
    end
end
