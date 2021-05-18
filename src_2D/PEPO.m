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
        testing
        visualise
        virtual_level_sizes_horiz
        virtual_level_sizes_vert
        PEPO_matrix
        boundary_matrix_x
        boundary_matrix_y
        boundary_vect
        bounds
        error_code
        copts
    end

    properties (Constant)
        numopts = struct('numbered', 1)
        cycleopts = struct('numbered', 1, ...
            'v_cyclic', 0, ...
            'h_cyclic', 1)
    end

    methods

        function [obj, err_code] = PEPO(model, opts, make_PEPO_handle)

            %parse opts
            p = inputParser;
            addParameter(p, 'order', 4) %max number of connected cells in chain. e.g. order 2 = 0--|--1--|--0
            addParameter(p, 'beta', 0.1)
            addParameter(p, 'testing', 0) %debug info
            addParameter(p, 'visualise', 0) %
            addParameter(p, 'inv_eps', 1e-12) %value for pseudo inverse
            addParameter(p, 'err_tol', 1e-13) %check whether constructed blocks are good enough
            addParameter(p, 'max_bond_dim', 20); %truncate bonds larger than this value
            addParameter(p, 'do_loops', 0); %parameters for 2D construction
            addParameter(p, 'loop_extension', 0);
            addParameter(p, 'double_extension', 0);
            addParameter(p, 'offset_loops', 0);
            addParameter(p, 'double_loop', 0);
            addParameter(p, 'complex', false); % complex or real PEPO construction

            parse(p, opts)

            order = p.Results.order;

            if mod(order, 2)
                max_index = (order - 1) / 2;
            else
                max_index = order / 2;
            end

            obj.copts = struct(...
                'L', max_index, ...
                'max_bond_dim', p.Results.max_bond_dim, ...
                'order', order, ...
                'inv_eps', p.Results.inv_eps, ...
                'complex', p.Results.complex, ...
                'err_tol', p.Results.err_tol, ...
                'do_loops', p.Results.do_loops, ...
                'loop_extension', p.Results.loop_extension, ...
                'double_extension', p.Results.double_extension, ...
                'offset_loops', p.Results.offset_loops, ...
                'double_loop', p.Results.double_loop ...
                );

            beta = p.Results.beta;

            obj.dim = model.d;
            obj.H_1_tensor = -beta * model.H_1_tensor;
            obj.H_2_tensor = -beta * model.H_2_tensor;

            ncells = 15; %make large enough for all possible simulations
            obj.PEPO_cell = cell(ncells, ncells, ncells, ncells);
            %obj.boundary_matrix_x = cell(max_index + 1, max_index + 2);
            %obj.boundary_matrix_y = cell(max_index + 1, max_index + 2);
            %obj.boundary_matrix_x{1, 1} = reshape(1, 1, 1);
            %obj.boundary_matrix_y{1, 1} = reshape(1, 1, 1);

            obj.virtual_level_sizes_horiz = 1;
            obj.virtual_level_sizes_vert = 1;

            obj.testing = p.Results.testing;
            obj.visualise = p.Results.visualise;

            %calculate ln of normailisation fac
            map = create_map(1:2, obj.numopts);
            [~, nf2] = H_exp(obj, map, 0, true);
            obj.nf = nf2;

            d = obj.dim;
            obj.PEPO_cell{1, 1, 1, 1} = reshape(eye(d) / exp(obj.nf), [d, d, 1, 1, 1, 1]);

            [obj, err_code] = make_PEPO_handle(obj);

            obj = cell2matrix(obj); %save matrix form

            obj.error_code = err_code;

            if err_code == 1
                warning("PEPO costruction failed, try lower beta");
            end

        end
    end
end
