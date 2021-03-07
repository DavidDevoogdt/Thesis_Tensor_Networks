function iorder = iorder(order)
    %function to inverse permutation vector
    n = numel(order);
    iorder = 1:n;
    iorder(order) = iorder;

end
