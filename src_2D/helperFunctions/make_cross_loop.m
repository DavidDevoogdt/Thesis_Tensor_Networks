function map = make_cross_loop(pat, m)
    %make array 1loop and 2 legs with given legth

    a = pat;
    a(m == 1) = 1;

    map = make_cross(a);

    if m(1) == 1
        x = 1;
    else
        x = size(map, 2);
    end

    if m(2) == 1
        y = 1;
    else
        y = size(map, 1);
    end

    map(y, x) = 1;

end
