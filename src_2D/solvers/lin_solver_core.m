function X = lin_solver_core(A_list, target)
    target_size_orig = size(target);

    num_A = numel(A_list);
    if numel(size(target)) ~= num_A +1
        target = reshape(target, [target_size_orig(1:num_A), prod(target_size_orig(num_A + 1:end))]);
    end

    R_dim = prod(size(target, 1:num_A));
    phys_dim = size(target, num_A + 1);

    %Q_cell = cell(num_A,1);
    B_cell = cell(1, num_A + 1);
    B_cell{1} = target;

    R_cell = cell(1, num_A);
    con_list = cell(1, num_A + 1);
    con_list{1} = [1:num_A, -(num_A + 1)];

    R_con_cell = cell(1, num_A);

    for i = 1:num_A
        [Q, R_cell{i}] = qr(A_list{i});
        con_list{i + 1} = [-i, i];
        R_con_cell{i} = [-i, -(num_A + i)];
        B_cell{i + 1} = Q';
    end

    B2 = ncon(B_cell, con_list);
    R = ncon(R_cell, R_con_cell);

    %B2 = ncon(  {Q1'*B1,Q2'*B2,Q3'*B3}, {  [-1,-4],[-2,-5],[-3,-6] }  );
    %R = ncon(  {R1,R2,R3}, {  [-1,-4],[-2,-5],[-3,-6] }  );

    R = reshape(R, [R_dim, R_dim]);
    B2 = reshape(B2, [R_dim, phys_dim]);

    dA = decomposition(R, 'triangular');

    X = dA \ B2;

    X = reshape(X, target_size_orig);

end
