function obj = cell2matrix(obj, sparse)
    %takes PEPO obj and puts all individual tensors in on big tensor

    if nargin < 2
        sparse = 0;
    end

    d = obj.dim;

    size_arr_horiz = obj.virtual_level_sizes_horiz;
    h_size = numel(size_arr_horiz) + 1;

    start_index_H = zeros(h_size, 1);
    start_index_H(1) = 1;
    ind = 1;

    for i = 2:h_size
        ind = ind + size_arr_horiz(i - 1);
        start_index_H(i) = ind;
    end

    totaldimensionH = start_index_H(end) - 1;

    function y = getH(i)
        y = start_index_H(i):start_index_H(i + 1) - 1;
    end

    size_arr_vert = obj.virtual_level_sizes_vert;
    v_size = numel(size_arr_vert) + 1;

    start_index_V = zeros(v_size, 1);
    start_index_V(1) = 1;
    ind = 1;

    for i = 2:v_size
        ind = ind + size_arr_vert(i - 1);
        start_index_V(i) = ind;
    end

    totaldimensionV = start_index_V(end) - 1;

    function y = getV(i)
        y = start_index_V(i):start_index_V(i + 1) - 1;
    end

    if sparse == 0
        T = zeros(d, d, totaldimensionH, totaldimensionV, totaldimensionH, totaldimensionV);
    else
        sparsem = ndSparse(sparse(d^2, totaldimensionH^2 * totaldimensionV^2));
        T = reshape(sparsem, [d, d, totaldimensionH, totaldimensionV, totaldimensionH, totaldimensionV]);
    end

    %move all existing tensors to matrix
    for i1 = 1:h_size - 1
        for i2 = 1:v_size - 1
            for i3 = 1:h_size - 1
                for i4 = 1:v_size - 1
                    %trace the spins for the environment
                    cell = obj.PEPO_cell{i1, i2, i3, i4};

                    if ~isempty(cell)
                        T(:, :, getH(i1), getV(i2), getH(i3), getV(i4)) = obj.PEPO_cell{i1, i2, i3, i4}; %ncon(  { obj.PEPO_cell{i1,i2,i3,i4} } , {[1,1,-1,-2,-3,-4]} );
                    end
                end
            end
        end
    end

    obj.PEPO_matrix = T;

    %K=reshape(T(:,1,1,:),[45,45]);

    %             obj.left = zeros(1, totaldimension);
    %             obj.left(1) = 1;
    %             obj.right = zeros(totaldimension, 1);
    %             obj.right(1) = 1;

end
