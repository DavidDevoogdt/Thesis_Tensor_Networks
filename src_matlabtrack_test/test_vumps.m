clear all; format long; doPath; close all

opts.charges='regular';
opts.dynamical='off';
opts.dyncharges=0;
opts.schmidtcut=1e-10;
opts.chimax=350;
opts.disp='iter';
opts.tolmax=1e-4;
opts.tolfactor=1e4;
opts.minit=1;
opts.dyniter=5;
opts.truncate=0;
opts.method='vumps';
opts.save=0;



opts.plot='on';
opts.maxit=1000;
opts.tolfixed=1e-8;
[O,M]=MpoIsing('none');

[A,G,lambda,ctr,error]=Vumps(O,40,[],opts);

%%

