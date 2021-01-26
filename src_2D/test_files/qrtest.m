n = 16;


A1 = rand(n,n);
A2 = rand(n,n);
A3 = rand(n,n);

A = ncon(  {A1,A2,A3}, {  [-1,-4],[-2,-5],[-3,-6] }  );

B1 = rand(n,1);
B2 = rand(n,1);
B3 = rand(n,1);

B = ncon(  {B1,B2,B3}, {  [-1,-4],[-2,-5],[-3,-6] }  );

A = reshape(A,[n^3,n^3]);
B = reshape(B,[n^3,1]);
%

tic
    dA = decomposition(A,'qr');%,'RankTolerance', 1e-12);
    X1 = dA\B;
toc

%
tic
    [U1,S1,V1] = svd(A1);
    [U2,S2,V2] = svd(A2);
    [U3,S3,V3] = svd(A3);

    U = ncon(  {U1,U2,U3}, {  [-1,-4],[-2,-5],[-3,-6] }  );
    V = ncon(  {V1,V2,V3}, {  [-1,-4],[-2,-5],[-3,-6] }  );
    S = ncon(  {S1,S2,S3}, {  [-1,-4],[-2,-5],[-3,-6] }  );

    U = reshape(U,[n^3,n^3]);
    S = reshape(S,[n^3,n^3]);
    V = reshape(V,[n^3,n^3]);

    S = diag(S);

    %[S_sort,idx] = sort(S,'descend');
    
    tol = max(size(A))*eps(norm(S));

    S( S<tol ) = 0;

    X2 =   V* diag(S.^(-1))*  U'*B;
    
toc

tic
    [Q1,R1] = qr(A1);
    [Q2,R2] = qr(A2);
    [Q3,R3] = qr(A3);

    B2 = ncon(  {Q1'*B1,Q2'*B2,Q3'*B3}, {  [-1,-4],[-2,-5],[-3,-6] }  );
    R = ncon(  {R1,R2,R3}, {  [-1,-4],[-2,-5],[-3,-6] }  );

    R = reshape(R,[n^3,n^3]);
    B2 = reshape(B2,[n^3,1]);
   
    dA = decomposition(R,'triangular');
    
    X3 = dA\B2;

toc

tic
    dA = decomposition(A,'cod');%,'RankTolerance', 1e-12);
    X4 = dA\B;
toc

n1 = svds( A*X1-B,1  );
n2 = svds( A*X2-B,1  );
n3 = svds( A*X3-B,1  );
n4 = svds( A*X4-B,1  );
fprintf("err1 %.4e\nerr2 %.4e\nerr3 %.4e\nerr4 %.4e\n",n1,n2,n3,n4);




