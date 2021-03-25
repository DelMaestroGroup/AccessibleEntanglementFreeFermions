# FreeFermionsOEE

A Julia  code for the modified correlation matrix method for spinless free fermions with a focus on operational entanglement entropy.

## Requirements
```julia
 pkg> add ArgParse
```
* https://github.com/carlobaldassi/ArgParse.jl

## Examples
```julia
 julia FFOEE_main.jl --l-step 1 --l-num 100 --l-min 1 200 100 --alpha 2 --out L200N100alpha2.dat --eigs --prob
 
 julia FFOEE_main.jl --precision 1000 --l-log --l-logstep 0.1 --l-min 5 --l-max 1000 1000 500 --alpha 2 --out L1000N500alpha2.dat
```
