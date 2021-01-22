function x_cell = split_x (x, x_sizes)

    num_x = size(x_sizes, 2);
    x_cell = cell(num_x, 1);

    curr = 0;

    for i = 1:num_x
        num_elem = prod(x_sizes{i});
        x_cell{i} = reshape(x(curr + 1:curr + num_elem), x_sizes{i});
        curr = curr + num_elem;
    end
end
