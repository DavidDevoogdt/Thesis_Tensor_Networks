function [obj,ln_prefact_out,err] = solve_lin_non_lin_and_assign(obj, map, patterns, ln_prefact, opts, assignfn, step)

    obj = fill_rand(obj, patterns, 1 / exp(obj.nf));

    if nargin >= 6
        obj = assignfn(obj);
    end

    p = inputParser;
    addParameter(p, 'maxit', 20)
    addParameter(p, 'solved', 1e-14)
    addParameter(p, 'minchange', 1 - 1e-4)
    addParameter(p, 'display', 0)
    parse(p, opts);

    maxit = p.Results.maxit;

    prev_err = Inf; %calculate_error(obj, map.num_map, obj.numopts);

    %precontract for speeed with matrix calculation and substract pattern
    %cells

    [all_con_cells, pat_cells] = get_valid_contractions(obj, map, struct('max_index', obj.current_max_index, 'pattern', {patterns}));

    con_cells = all_con_cells(pat_cells);
    pat_cells = pat_cells(pat_cells == 1);

    obj = cell2matrix(obj);
    [target, ln_prefact_out] = H_exp(obj, map, ln_prefact, true);
    target = target - contract_network(obj, map, struct('matrix', 1, 'lnprefact', ln_prefact_out));
    target_site = reshape(permute(target, site_ordering_permute(map.N)), dimension_vector(obj.dim^2, map.N));

    target_site = -contract_con_cells(obj, map, ln_prefact_out, -target_site, con_cells);

    %if obj.testing == 1
    %    fprintf("starting solver")
    %end

    % if nargin >= 6
    %     pcells = find(pat_cells);

    %     matching = zeros(size(pcells));
    %     for i = 1:numel(pcells)
    %         ccell = all_con_cells{pcells(i)}{1};

    %         match = zeros(size(ccell));
    %         for j = 1:numel(ccell)
    %             for k = 1:numel(patfilt)
    %                 if sum(patfilt{k} == ccell{j}) == 4
    %                     match(j) = 1;
    %                     break;
    %                 end
    %             end
    %         end
    %         if match == ones(size(match))
    %             matching(i) = 1;
    %         end
    %     end

    %     pat_cells(pcells(~matching)) = 0;

    % end

    for j = 1:maxit

        for i = 1:numel(patterns)
            if nargin >= 6
                leg = patterns{1};
                prev_A = obj.PEPO_cell{leg(1) + 1, leg(2) + 1, leg(3) + 1, leg(4) + 1};
            end
            
            obj = solve_lin_and_assign(obj, map, patterns(i), ln_prefact_out, struct('target_site',target_site ,'all_con_cells',{con_cells} ,'pat_cells',pat_cells));

        end

        if nargin < 6
            err_vect = contract_con_cells(obj, map, ln_prefact_out, target_site, con_cells);
            mul_factor = exp(ln_prefact_out - obj.nf);
            err = sum(svds(reshape(err_vect, [obj.dim^(map.N), obj.dim^(map.N)]), 10).^2).^0.5 * mul_factor;

        else

            leg = patterns{1};
            new = obj.PEPO_cell{leg(1) + 1, leg(2) + 1, leg(3) + 1, leg(4) + 1};
            old = prev_A;

            dist = norm(reshape(new - old, [], 1));

            ss = step / (dist);

            obj.PEPO_cell{leg(1) + 1, leg(2) + 1, leg(3) + 1, leg(4) + 1} = old + step * (new - old);

            %obj = rescale_PEPO_pattern(obj, full_pat,1.0001,0);

            obj = assignfn(obj);

            err = calculate_error(obj, map, [], 1);
        end

        if p.Results.display == 1
            fprintf("%3d: %.4e\n", j, err);
        end

        if err < p.Results.solved
            if p.Results.display == 1
                fprintf('under threshold')
            end
            break;
        end

        if err / prev_err > p.Results.minchange
            if nargin >= 6
                step = step * 0.9;

                if step < 0.05
                    fprintf("not solved, res err %.4e\n", err)
                    break;
                end

            else
                fprintf("not solved, res err %.4e\n", err)
                break;
            end
        end

        prev_err = err;
    end

end

function F = simil(x, A, B)
    F = ncon({inv(x), A, x}, {[-3, 1], [-1, -2, 1, 2], [2, -4]}) - B;
end
