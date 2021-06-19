# Thesis_Tensor_Networks
This git repo shows some of the work done for my master thesis at UGent quantum group http://mathphy.ugent.be/wp/quantum/

Report: https://github.com/DavidDevoogdt/Thesis_Tensor_Networks/blob/master/verslag/main.pdf

Presentation: https://github.com/DavidDevoogdt/Thesis_Tensor_Networks/blob/master/verslag/Eindpresentatie/main.pdf



## How to use

start Matlab in root folder of this repo, execute

```
doPath.m
test 
test_2D(struct('do_loops',1,'testing',1),1,'t_ising') 
Ising2D_par(8, 2.5, 'g', struct('testing',1,'unit_cell',1,'par',0,'order',5,'do_loops',1 ));
```
