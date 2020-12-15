
target_mps=  rand(4,4,4,4);
dim=6;

opts.tol = 1e-13;
opts.maxit = 5000;
opts.print_level = 1;

[elems,err] = tensor_ring(target_mps, dim,opts );

