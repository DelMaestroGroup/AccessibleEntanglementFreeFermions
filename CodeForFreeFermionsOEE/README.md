# FreeFermionsOEE

A Julia  code for the modified correlation matrix method for spinless free fermions with a focus on operational entanglement entropy.

## Requirements

* [ArgParse](https://github.com/carlobaldassi/ArgParse.jl) (`Pkg.add("ArgParse")`)

## Examples

* `julia tV_main.jl --help`
* `julia FFOEE_main.jl  --pbc  --l-step 1 --l-num 100  --l-min 1  200 100  --alpha 2 --out output.dat`
* `julia FFOEE_main.jl  --pbc  --precision 1000 --l-log --l-logstep 0.1 --l-min 5 --l-max 1000  1000 500  --alpha 2 --out output.dat`