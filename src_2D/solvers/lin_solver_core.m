function X = lin_solver_core(A_list, target, inv_eps)

    %R2_dims = []

    target_size_orig = size(target);

    num_A = numel(A_list);
    if numel(size(target)) ~= num_A +1
        target = reshape(target, [target_size_orig(1:num_A), prod(target_size_orig(num_A + 1:end))]);
    end

    type = 0;

    if type == 1%qr

        R_dim = prod(size(target, 1:num_A));
        phys_dim = size(target, num_A + 1);

        %Q_cell = cell(num_A,1);
        B_cell = cell(1, num_A + 1);
        B_cell{1} = target;

        R_cell_dims = zeros(1, num_A);

        R_cell = cell(1, num_A);
        con_list = cell(1, num_A + 1);
        con_list{1} = [1:num_A, -(num_A + 1)];

        R_con_cell = cell(1, num_A);

        for i = 1:num_A
            [Q, R_cell{i}] = qr(A_list{i});

            R_cell_dims(i) = size(A_list{i}, 2);

            con_list{i + 1} = [-i, i];
            R_con_cell{i} = [-i, -(num_A + i)];
            B_cell{i + 1} = Q';
        end

        B2 = ncon(B_cell, con_list);
        R = ncon(R_cell, R_con_cell);

        %B2 = ncon(  {Q1'*B1,Q2'*B2,Q3'*B3}, {  [-1,-4],[-2,-5],[-3,-6] }  );
        %R = ncon(  {R1,R2,R3}, {  [-1,-4],[-2,-5],[-3,-6] }  );

        R = reshape(R, R_dim, []);
        R_dim2 = size(R, 2);

        B2 = reshape(B2, [R_dim, phys_dim]);

        %B2(B2<1e-15)=0;

        t = diag(R);
        tol = max(size(R)) * eps(norm(t));
        tol = max(tol, inv_eps);

        mask = abs(diag(R)) < tol;
        if sum(mask) ~= 0
            sum(mask);
        end

        X = zeros(R_dim2, phys_dim);
        R2 = decomposition(R(~mask, ~mask), 'triangular');
        X(~mask, :) = R2 \ B2(~mask, :);

        %fprintf('%.4e ', max(max(abs(X))));

        %X = reshape(X, sz);

        X = reshape(X, [R_cell_dims, phys_dim]);

    else %svd

        S_dim = prod(size(target, 1:num_A));
        phys_dim = size(target, num_A + 1);

        %Q_cell = cell(num_A,1);
        U_cell = cell(1, num_A + 1);
        U_cell{1} = target;

        V_cell = cell(1, num_A);
        V_cell{1} = target;

        S_cell_dims = zeros(1, num_A);
        S_cell = cell(1, num_A);

        con_list = cell(1, num_A + 1);
        con_list{1} = [1:num_A, -(num_A + 1)];

        S_con_cell = cell(1, num_A);

        for i = 1:num_A
            [U, S, V_cell{i}] = svd(A_list{i});

            S_cell{i} = ndSparse(S);

            S_cell_dims(i) = size(A_list{i}, 2);

            con_list{i + 1} = [-i, i];
            S_con_cell{i} = [-i, -(num_A + i)];
            U_cell{i + 1} = U';
        end

        RHS = ncon(U_cell, con_list);
        S = ncon(S_cell, S_con_cell);

        %B2 = ncon(  {Q1'*B1,Q2'*B2,Q3'*B3}, {  [-1,-4],[-2,-5],[-3,-6] }  );
        %R = ncon(  {R1,R2,R3}, {  [-1,-4],[-2,-5],[-3,-6] }  );

        S = reshape(S, S_dim, []);
        RHS = reshape(RHS, [S_dim, phys_dim]);

        dS = decomposition(sparse(S), 'qr', 'RankTolerance', inv_eps, 'CheckCondition', false);
        Y = dS \ RHS;

        Y = reshape(Y, [S_cell_dims, phys_dim]);

        X = ncon([Y, V_cell], con_list);

    end

end
