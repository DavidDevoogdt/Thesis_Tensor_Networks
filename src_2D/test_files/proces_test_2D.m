% function test
%     %test_H_exp
%     test_PEPO
%     %test3
% end

%function test_PEPO

function test

    %model = 't_ising';
    model = 'heis';

    switch model
        case 't_ising'

            filenames = {
                '2D_05-20-21_18-43-17';%'2D_05-20-21_16-02-11';
                         %'2D_05-20-21_16-01-39';
                         '2D_05-20-21_16-42-07';
                         %'2D_05-20-21_16-42-11';
                         %'2D_05-20-21_17-25-55';
                         '2D_05-20-21_19-37-08';
                         '2D_05-20-21_20-27-31';
                         '2D_05-21-21_08-59-27'
                         '2D_05-21-21_13-54-44';
                         };

            labels = {
                %'$\epsiolon^1$ no loops, $\sigma_0=10^{-12}$';
                     '$\epsilon^1$ no loops';
                     '$\epsilon^1$ plaquette';
                     %'no loops, $\sigma_0=10^{-13}$';
                     %'$\epsiolon^1$plaquette, $\sigma_0=10^{-10}$';
                     '$\epsilon^1$ extensions';
                     '$\epsilon^2$ no loops';
                     '$\epsilon^2$ plaquette';
                     '$\epsilon^2$ extensions';
                     };
        case 'heis'
           filenames = {
               '2D_05-24-21_08-55-47';
               '2D_05-24-21_09-30-36';
               '2D_05-24-21_10-15-02';
               '2D_05-24-21_11-43-54';
               '2D_05-24-21_17-58-43';
               '2D_05-25-21_10-47-32';
               };
                     
                     
            labels = {
                     '$\epsilon^1$ no loops';
                     '$\epsilon^1$ plaquette';
                     '$\epsilon^1$ extensions';
                     '$\epsilon^2$ no loops';
                     '$\epsilon^2$ plaquette';
                     '$\epsilon^2$ extensions';
                     };

        otherwise
            error(unknow)
    end

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
