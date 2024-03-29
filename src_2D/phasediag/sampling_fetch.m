function [data, template] = sampling_fetch(file_name, opts)
    %helper function to read save files with a given name
    if ~isfield(opts, 'save_vars')
        opts.save_vars = 0;
    end

    do_call_back = isfield(opts, 'call_back_fn');

    fold = mfilename('fullpath');
    pathparts = strsplit(fold, '/');
    pathparts = [pathparts(1:end - 3), 'IsingMatFiles'];
    fold2 = strjoin(pathparts, '/');

    name = sprintf("%s.mat", file_name);
    folder = sprintf("%s/%s/", fold2, file_name);

    if ~exist(folder, 'dir')
        error("invalid files")
    end

    myFiles = dir(fullfile(folder, 'template.mat'));
    if numel(myFiles) ~= 1
        error('template not found')
    else
        baseFileName = myFiles(1).name;
        fullFileName = fullfile(folder, baseFileName);
        load(fullFileName, 'template');
        template.found = 1;
    end

    myFiles = dir(fullfile(folder, 'results_*.mat'));
    num_files = numel(myFiles);

    data_points = cell(num_files, 1);
    if opts.save_vars == 1
        save_data_points = cell(num_files, 1);
    end

    mergestructs = @(x, y) cell2struct([struct2cell(x); struct2cell(y)], [fieldnames(x); fieldnames(y)]);

    parfor k = 1:num_files

        baseFileName = myFiles(k).name;
        fullFileName = fullfile(folder, baseFileName);
        R = load(fullFileName, 'results');
        data_points{k} = R.results;
            
        S = [];
        if opts.save_vars == 1
            %baseFileName = myFiles(k).name;
            baseFileName = strrep(baseFileName, 'results', 'save_vars');
            fullFileName = fullfile(folder, baseFileName);

            S = load(fullFileName, 'save_vars');

            save_data_points{k} = S.save_vars;
        end

       
        
        if do_call_back
            
            opts.call_back_fn(R.results, S.save_vars, baseFileName);
        end
    end
    
    if nargout ~= 0

        fields = fieldnames(data_points{1});

        reordered_data = struct();

        %valid = cell2mat(cellfun(@(x) x.err.', data_points, 'UniformOutput', false)) < template.vumps_opts.tolfixed;
        %data_points = data_points(valid);
        
        
        for i = 1:numel(fields)
            fn = fields{i};

            if fn== "eps_i"
                val = cellfun(@(x) x.(fn)(1:6).', data_points, 'UniformOutput', false);
            else
                val = cellfun(@(x) x.(fn).', data_points, 'UniformOutput', false);
            end
            
            
            reordered_data.(fn) = cell2mat(val);
        end

        if opts.save_vars == 1
            fields2 = fieldnames(save_data_points{1});
            for i = 1:numel(fields2)
                fn = fields2{i};
                val = cellfun(@(x) x.(fn), save_data_points, 'UniformOutput', false);
                reordered_data.(fn) = val;
            end
        end

        %add template stuff

        data = mergestructs(reordered_data, template);

        if opts.save_vars == 1
            data.fields = [fields; fields2];
        else
            data.fields = fields;
        end

        %complete incomplete files

        if ~isfield(data, 'free_var')
            data.free_var = 'T';
            data.fixed_var = 'g';
            data.fixed_val = data.model_params.g;
        end
    end

end
