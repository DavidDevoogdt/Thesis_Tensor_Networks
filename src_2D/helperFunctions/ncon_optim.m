%uses netcon to generate the optimal sequence and returns results from
%ncon;

function res = ncon_optim(tensor_list, leg_list)

    pos_list = [];
    neg_list = [];
    counter = 1;

    for i = 1:size(leg_list, 2)

        list = leg_list{i}; %set size trailing elements to 1

        for j = 1:size(list, 2)

            if list(j) > 0
                pos_list(list(j), :) = [list(j), size(tensor_list{i}, j)];
            else
                neg_list(counter, :) = [list(j), size(tensor_list{i}, j)];
                counter = counter + 1;
            end

        end

    end

    cost_list = cat(1, pos_list, neg_list);

    [seq] = netcon(leg_list, 0, 1, 1, 1, cost_list);

    res = ncon(tensor_list, leg_list, seq);
end
