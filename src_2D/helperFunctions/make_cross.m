function map = make_cross(arr)
    mx = arr(1) + 1;
    my = arr(2) + 1;

    x = arr(1) + arr(3) + 1;
    y = arr(2) + arr(4) + 1;

    map = zeros([y, x]);

    map(my, mx) = 1;

    map(my, 1:arr(1)) = 1;
    map(my, mx:mx + arr(3)) = 1;
    map(1:arr(2), mx) = 1;
    map(my:my + arr(4), mx) = 1;
end
