function tensors = fetch_PEPO_cells(obj, map, legs, ln_prefactor, patterns, xs)

    if nargin < 4
        ln_prefactor = obj.nf;
    end

    mult_fact = exp(ln_prefactor - obj.nf);

    if nargin < 5
        patterns = [];
    end

    num_patterns = size(patterns, 2);

    tensors = cell(1, map.N);

    for n = 1:map.N

        leg = legs{n};

        matched_pattern = 0;

        for ii = 1:num_patterns

            if same_pattern(leg, patterns{ii}) == 1
                tensorsn = xs{ii};
                matched_pattern = 1;
                break;
            end

        end

        if matched_pattern == 0

            if map.is_x_border(n)
                tensorsn = obj.boundary_matrix_x{leg(1) + 1, leg(2) + 1};
            elseif map.is_y_border(n)
                tensorsn = obj.boundary_matrix_y{leg(1) + 1, leg(2) + 1};
            else
                tensorsn = obj.PEPO_cell{leg(1) + 1, leg(2) + 1, leg(3) + 1, leg(4) + 1};
            end

        end

        tensors{n} = tensorsn / mult_fact;
    end

end