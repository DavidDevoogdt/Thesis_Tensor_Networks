% function test
%     %test_H_exp
%     test_PEPO
%     %test3
% end

%function test_PEPO

function test

    filenames = {'2D_05-20-21_18-43-17';%'2D_05-20-21_16-02-11';
                 '2D_05-20-21_16-01-39';
                 '2D_05-20-21_16-42-07';
                 %'2D_05-20-21_16-42-11';
                 '2D_05-20-21_17-25-55';
                 '2D_05-20-21_19-37-08';};
             
    labels = {'no loops, $\sigma_0=10^{-12}$';
             'no loops, $\sigma_0=10^{-10}$';
             'plaquette, $\sigma_0=10^{-12}$';
             %'no loops, $\sigma_0=10^{-13}$';
             'plaquette, $\sigma_0=10^{-10}$';
             'extensions'};

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');

    pathparts = [pathparts(1:end - 3), 'test_2D_files'];
    fold2 = strjoin(pathparts, '/');
    
    
    figure();
    x_width = 15;
    y_width = 10;
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, x_width, y_width], 'PaperSize', [x_width, y_width],'DefaultAxesFontSize', 11)

    for i =1:numel(filenames)

        filename = sprintf("%s/%s.mat", fold2, filenames{i});
        fprintf("%s \n", filename)
        load(filename, 'time', 'simul', 'beta_arr', 'err_arr', 'num_map', 'map_opts', 'density_site');


        loglog(beta_arr, err_arr);
      
        hold on

        if i==1
            title(simul.title);
            xlabel("$\beta \cdot J $", "Interpreter", "Latex");
            ylabel("err");
            ylim([1e-3, 1e2])
            ylim([0, 10])
        end

        figure(gcf)
    
    end
    
    
    legend(labels, 'Location', 'northwest',"Interpreter", "Latex",'Fontsize',12)
    
    hold off
    

end
