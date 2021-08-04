function doPath()

    % add all subfolders to path
    % author: LV

    warning('off', 'MATLAB:mpath:nameNonexistentOrNotADirectory')

    addpath(genpath(pwd));

end
