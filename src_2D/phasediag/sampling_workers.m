function sampling_workers(aim_dx, aim_dy, nsample, template, npoints, first, par)

    if nargin < 7
        par = 1;
    end

    x_min = template.x_bounds(1);
    x_max = template.x_bounds(2);

    if par == 0

        save_vars = [];
        save_vars.fname = sprintf('%2d:%2d', 1, 1);
        v = sampling_core(save_vars, template, x_min);
        fprintf(v);

        save_vars = [];
        save_vars.fname = sprintf('%2d:%2d', template.iter, 2);
        v = sampling_core(save_vars, template, x_max);
        fprintf(v);

        counter = 3;

        while counter < npoints

            data = sampling_fetch(template.name, struct);
            data = sampling_filter(data, struct);
            %determine new point
            x0 = get_new_point(data.(template.free_var), data.m, [], ratio);

            save_vars = [];
            save_vars.fname = sprintf('%2d:%2d', template.iter, counter);

            v = sampling_core(save_vars, template, x0, 1);
            fprintf(v);

            counter = counter + 1;
        end
    else
        f(1:nsample) = parallel.FevalFuture;

        counter = nsample + 1;
        ratio = aim_dx / aim_dy;

        if first == 1
            x0 = (x_max - x_min) / (nsample - 1) * (0:nsample - 1) + x_min;

            for j = 1:nsample
                save_vars = [];
                save_vars.fname = sprintf('%2d:%2d', template.iter, j);

                f(j) = parfeval(@ () sampling_core(save_vars, template, x0(j)), 1);
            end
        else
            x0 = [];
            for j = 1:nsample
                x0 = [x0, 0];
                [f, x0] = add_new_point(template, j, f, ratio, x0, j);
            end
        end

        while counter < npoints
            [idx, v] = fetchNext(f);
            fprintf(v);

            if counter < npoints - nsample
                [f, x0] = add_new_point(template, idx, f, ratio, x0, counter);
                counter = counter + 1;
            end
        end
    end
end

function [f, x0] = add_new_point(template, idx, f, ratio, x0, number)
    data = sampling_fetch(template.name, struct);
    data = sampling_filter(data, struct);

    %determine new point
    z0 = x0;
    z0(idx) = [];

    x0(idx) = get_new_point(data.(template.free_var), data.m, z0, ratio);

    save_vars = [];
    save_vars.fname = sprintf('%2d:%2d', template.iter, number);

    f(idx) = parfeval(@ () sampling_core(save_vars, template, x0(idx)), 1);
end

function x = get_new_point(x, y, x0, ratio)

    x_arr = [x.', x0];

    if numel(x) == 1
        y_arr = [y, x0 * 0 + y];
    else
        y_arr = interp1(x, y, x_arr, 'linear', 'extrap');
    end

    [x_arr, I] = sort(x_arr);
    y_arr = y_arr(I);

    dx_arr = diff(x_arr);
    dy_arr = diff(y_arr);

    ds = (dx_arr.^2 + (ratio .* dy_arr).^2).^(0.5);
    [~, idx] = max(ds);
    x = x_arr(idx) + dx_arr(idx) / 2;
end

function x0 = get_new_points(data, aim_dx, aim_dy, template)

    error('check this function')

    data = sampling_fetch(template.name, struct);
    data = sampling_filter(data, struct);

    x_arr = data.(template.free_var);
    y_arr = data.m;

    dx_arr = diff(x_arr);
    dy_arr = diff(y_arr);

    ds = (dx_arr.^2 + ((aim_dx / aim_dy) .* dy_arr).^2).^(0.5);

    new_x = zeros(size(ds));

    for ii = 1:nsample
        [~, idx] = max(ds ./ (new_x + 1));
        new_x(idx) = new_x(idx) + 1;
    end

    x0 = [];

    K = find(new_x ~= 0);
    for ii = 1:numel(K)
        ind = K(ii);
        nnx = new_x(ind);
        nx = x_arr(ind);
        dx = dx_arr(ind) / (nnx + 1);

        x0 = [x0, nx + dx * (1:nnx)];
    end
end
