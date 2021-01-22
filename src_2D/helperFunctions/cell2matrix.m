function obj = cell2matrix(obj)

    d = obj.dim;

    size_arr_horiz = obj.virtual_level_sizes_horiz;

    start_index_H = zeros(obj.max_index + 2, 1);
    start_index_H(1) = 1;
    ind = 1;

    for i = 2:obj.max_index + 2
        ind = ind + size_arr_horiz(i - 1);
        start_index_H(i) = ind;
    end

    totaldimensionH = start_index_H(end) - 1;

    function y = getH(i)
        y = start_index_H(i):start_index_H(i + 1) - 1;
    end

    size_arr_vert = obj.virtual_level_sizes_vert;

    start_index_V = zeros(obj.max_index + 2, 1);
    start_index_V(1) = 1;
    ind = 1;

    for i = 2:obj.max_index + 2
        ind = ind + size_arr_vert(i - 1);
        start_index_V(i) = ind;
    end

    totaldimensionV = start_index_V(end) - 1;

    function y = getV(i)
        y = start_index_V(i):start_index_V(i + 1) - 1;
    end

    T = zeros(d, d, totaldimensionH, totaldimensionV, totaldimensionH, totaldimensionV);

    %sparsem = ndSparse(sparse(d^2, totaldimensionH^2*totaldimensionV^2));
    %T = reshape(sparsem, [d,d,totaldimensionH,totaldimensionV,totaldimensionH,totaldimensionV]);

    %move all existing tensors to matrix
    for i1 = 1:obj.max_index + 1

        for i2 = 1:obj.max_index + 1

            for i3 = 1:obj.max_index + 1

                for i4 = 1:obj.max_index + 1
                    %trace the spins for the environment
                    cell = obj.PEPO_cell{i1, i2, i3, i4};

                    if length(cell) ~= 0
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