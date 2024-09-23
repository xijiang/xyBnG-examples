# Package xyBnG examples

## Introduction

This repository contains examples of the package `xyBnG`. The package `xyBnG` is a Julia package for simulations of complicate breeding schemes. The package xyBnG is available at https://github.com/xijiang/xyBnG.jl.

I haven't register it as a standard Julia package yet. This will be done when I feel comfortable with the package. You can add it by the URL of the package repository. 

This repo, i.e., `xyBnG_examples`, contains examples of using the package `xyBnG`. One need to install the package `xyBnG` before running these examples.

The files here come and go. 'Good' examples will be moved to the package repository.

## `pgsnp.jl`

This script read results from `pgsnp.cpp` as founder population and simulate a breeding scheme with the package `xyBnG`. 

### `pgsnp.cpp`

```bash
g++ -Wall -O3 -std=c++11 -pthread -o pgsnp pgsnp.cpp
ulimit -s 600000 # increase stack size when necessary
./pgsnp 200 2000 1.5 1 > snp.txt
```

This program was written in 2010. It uses too much stack memory. The `ulimit -s 600000` is to increase the stack size.

### `pgsnp.jl`

Warning: This script calls `pgsnp`. You have to make sure that the stack size if large enough. 

```bash
ulimit -s 600000 # increase stack size when necessary
julia
```

You are also recomended to have `Revise` installed. Insert the following code in `~/.julia/config/startup.jl`:

```julia
using Revise
```

Above will load `Revise` automatically when you start Julia. In the future, if you `includet("your-script.jl")`, the script will be revised automatically when you save it. You don't have to restart Julia, or `includet` your scripts again.

Run the following code in Julia to us the script `pgsnp.jl`:

```julia
includet("path-to/pgsnp.jl")
pgsnp()
```
