% function test
%     %test_H_exp
%     test_PEPO
%     %test3
% end

%function test_PEPO

function test

    filename = "2D_03-21-21_10-33-49";

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');

    pathparts = [pathparts(1:end - 3), 'test_2D_files'];
    fold2 = strjoin(pathparts, '/');

    filename = sprintf("%s/%s.mat", fold2, filename);

    fprintf("%s \n", filename)

    load(filename)

    %%

    figure();
    x_width = 15;
    y_width = 10;
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, x_width, y_width], 'PaperSize', [x_width, y_width])

    loglog(beta_arr, err_arr);
    %plot(T,err_arr );

    %hold on

    title(sprintf("2D transverse Ising, g=%.4f ", g));
    xlabel("$\frac{J}{k T}$", "Interpreter", "Latex");
    ylabel("err");
    ylim([0, 10])

    figure(gcf)

end
