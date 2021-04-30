function [map,pat] = make_cross_loop(arr,li)
    mx = arr(1) + 1;
    my = arr(2) + 1;

    x = arr(1)  + 2;
    y = arr(2)  + 2;

    map = zeros([y, x]);

    map(my:my+1, mx:mx+1) = 1;

    map(my, 1:arr(1)) = 1;
    map(1:arr(2), mx) = 1;
    
    pat = [arr(1),arr(2),li(1),li(2)];
end
