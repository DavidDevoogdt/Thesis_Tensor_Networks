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
        nf%normalisation factor
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
        order
    end

    methods

        function obj = PEPO(d, H_1_tensor, H_2_tensor, order, make_PEPO_handle, opts)
            numopts.numbered = 1;
            obj.numopts = numopts;

            cycleopts.numbered = 1;
            cycleopts.v_cyclic = 0;
            cycleopts.h_cyclic = 1;
            obj.cycleopts = cycleopts;

            obj.dim = d;
            obj.H_1_tensor = H_1_tensor;
            obj.H_2_tensor = H_2_tensor;

            if mod(order, 2)
                max_index = (order - 1) / 2;
            else
                max_index = order / 2;
            end
            obj.max_index = max_index;
            obj.order = order;
            obj.cycle_index = Inf;

            obj.PEPO_cell = cell(max_index + 1, max_index + 1, max_index + 1, max_index + 1);
            obj.boundary_matrix_x = cell(max_index + 1, max_index + 1);
            obj.boundary_matrix_y = cell(max_index + 1, max_index + 1);
            obj.boundary_matrix_x{1, 1} = reshape(1, 1, 1);
            obj.boundary_matrix_y{1, 1} = reshape(1, 1, 1);
            obj.virtual_level_sizes_horiz = 1;
            obj.virtual_level_sizes_vert = 1;

            %parse opts
            p = inputParser;
            addParameter(p, 'testing', 0)
            addParameter(p, 'visualise', 0)
            parse(p, opts)

            obj.testing = p.Results.testing;
            obj.visualise = p.Results.visualise;

            %calculate ln of normailisation fac
            map = create_map(1:2, numopts);
            [~, nf2] = H_exp(obj, map, 0, true);
            obj.nf = nf2;

            obj.PEPO_cell{1, 1, 1, 1} = reshape(eye(d) / exp(obj.nf), [d, d, 1, 1, 1, 1]);

            %non generic PEPO code
            obj = make_PEPO_handle(obj);

            obj = cell2matrix(obj); %save matrix form

        end
    end
end
