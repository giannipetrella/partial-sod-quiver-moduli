# Verify the conjecture for a given case

This repository is part of the paper "Partial semiorthogonal decomposition
for quiver moduli" of Gianni Petrella.

The code in this repository allows to verify whether Questions A, B and C as described
in the paper hold for any given datum of a quiver `Q` and any dimension vector `d`.

The script `verification.jl` can reproduce all the results for which it is cited in the paper.

To use this code, follow the instructions below.


## Instructions

1. You need a Julia installation.
Install Julia from [the Julia download page](https://julialang.org).

2. Run the Julia REPL in the folder where you wish to save `verification.jl` and
install `QuiverTools` by running

```julia
using Pkg; Pkg.add(url="https://github.com/QuiverTools/QuiverTools.jl.git")
```

3. Quit the Julia REPL and in your terminal run

```sh
julia verification.jl
```

To verify whether any of the questions of the paper holds for your favourite quiver `Q`,
dimension vector `d` and choice of linearisation `a`,
comment out the examples from the paper in `verification.jl` and replace them with yours.
If you don't want to guess a linearisation `a`, `QuiverTools` can find one: use

```julia
a = QuiverTools.extended_gcd(d)[2]
```

This work was supported by the Luxembourg National Research Fund (FNR–17953441).
