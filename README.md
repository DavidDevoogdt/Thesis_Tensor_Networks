# Thesis_Tensor_Networks
This git repo shows some of the work done for my master thesis at UGent quantum group http://mathphy.ugent.be/wp/quantum/

unfinished report: https://github.com/DavidDevoogdt/Thesis_Tensor_Networks/blob/master/verslag/main.pdf


How To use:

start matlab in root folder of this repo, excute

```
doPath.m
test 
test_2D(struct('do_loops',1,'testing',1),1,'t_ising') 
Ising2D_par(8, 2.5, 'g', struct('testing',1,'unit_cell',1,'par',0) , struct('order',5 ,'do_loops',1 ));
```
