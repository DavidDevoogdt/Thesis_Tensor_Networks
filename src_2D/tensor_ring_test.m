% mapopts.numbered=0;
% mapopts.v_cyclic=1;
% mapopts.h_cyclic=1;
% 
% map = PEPO.create_map([0,0,1,1;
%                        1,1,1,1;
%                        0,0,1,1;],mapopts); %cyclic ring
% 

% mapopts.numbered=0;
% mapopts.v_cyclic=0;
% mapopts.h_cyclic=1;
%map = PEPO.create_map([1,1,1,1],mapopts); %cyclic ring


%fit random 1 loop

mapopts.numbered=1;
mapopts.v_cyclic=0;
mapopts.h_cyclic=0;

map = PEPO.create_map([1,2;
                       3,4],mapopts); %loop

physdim = 2;

target_mps=  rand(physdim^2,physdim^2,physdim^2,physdim^2);
bond_dim=6;

opts.tol = 1e-13;
opts.maxit = 5000;
opts.print_level = 1;


elem_list = cell(1,4);

bond_dim_h=9;
bond_dim_v=9;

elem_list{1} = rand(  physdim,physdim, 1,1,bond_dim_h,bond_dim_v);
elem_list{2} = rand(  physdim,physdim, bond_dim_h,1,1,bond_dim_v);
elem_list{3} = rand(  physdim,physdim, 1,bond_dim_v,bond_dim_h,1);
elem_list{4} = rand(  physdim,physdim, bond_dim_h,bond_dim_v,1,1);


[elems,err] = tensor_ring(elem_list,map,target_mps, opts);

%firt double lopp

target_mps_2=  rand(physdim^2,physdim^2,physdim^2,physdim^2,physdim^2,physdim^2);

opts2.tol = 1e-13;
opts2.maxit = 5000;
opts2.print_level = 1;
opts2.optim = [5,6];

elem_list_2 = cell(1,6);




elem_list_2(1:4) = elems(1:4);
elem_list_2{5} = rand(  physdim,physdim, bond_dim_h,1,bond_dim_h,bond_dim_v);
elem_list_2{6} = rand(  physdim,physdim, bond_dim_h,bond_dim_v,bond_dim_h,1);

map = PEPO.create_map([1,5,2;
                       3,6,4],mapopts); 
[elems2,err2] = tensor_ring(elem_list_2,map,target_mps_2, opts2);               
                   
                   



