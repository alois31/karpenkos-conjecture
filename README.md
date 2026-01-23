# Counterexamples to Karpenko's conjecture for Spin groups

We aim to construct counterexamples to Karpenko's conjecture for Spin groups, by finding torsion in the module of irrational elements of the truncated Brown-Peterson cohomology of the maximal oriented Grassmannian.
The required computation are performed using a computer with the help of the present set of programs.

## Prerequisites

The software has only been tested on Linux operating systems.
Other platforms are expected to work as well as long as all dependencies are available.

* The [Rust toolchain](https://rust-lang.org/tools/install/).
  Most testing has been performed with version 1.86.0, but in principle all versions starting from 1.82.0 should work.

* If wishing to perform exact computations (see below), the [Singular](https://www.singular.uni-kl.de/) computer algebra system for polynomial computations.
  At least version 4.4.1 with additional, at the time of writing unreleased, patches ([1](https://github.com/Singular/Singular/commit/926ff5410741911f311f3cef9ddfb4fd7c789bc3), [2](https://github.com/Singular/Singular/commit/0e16e23693fa02fbdc5e5208489b18a15fe30318), [3](https://github.com/Singular/Singular/commit/ac5186c27ddaa91a04cdd7a294fb91bdfc0ba281)) is required.

  **Warning:** Do not attempt to use unpatched version 4.4.1, as it contains a bug (fixed by the patches) leading to miscomputations.
  Older versions contain at least one additional bug also leading to miscomputations, so will not work even with the patches applied.

If the [Lix](https://lix.systems/install/) or [Nix](https://nixos.org/download/) package manager is available, the easiest way to obtain the required software is to enter the shell environment using `nix-shell --arg development false`.

## Exact computations

For n≤7, the Brown-Peterson cohomology of the maximal orthogonal Grassmannian can be computed exactly.
The computations will be quite slow, particularly when n is large (around 1 hour for n=7).
In the code example below, the shell variable `n` needs to be set to the desired number.

    cargo run -p generate-exact --release -- $n
    Singular -q results/ogr$n/stage1.sing > results/ogr$n/stage2.sing
    Singular -q results/ogr$n/stage3.sing

## Approximate computations

Larger n are only handled approximately up to O(v^5).
Theoretically all n≤22 are supported, but n≤13 is recommended due to excessive running time, particularly for n≥17 where things would get interesting again.
Again, computation will take some time, especially for large n (around half a day for n=13).
The shell variable `n` needs to be set to the desired number again.

    cargo run -p approximate --release -- $n

## Pre-computed results

Due to the long running time of the programs, for convenience pre-computed results are available in the separate [results](https://codeberg.org/alois3264/karpenkos-conjecture/src/branch/results) branch.

## Citation

The paper is not published yet, final citation information will be inserted later.

    @misc{karpenkos-conjecture-counterexamples,
        author = {Petrov, V. and Wohlschlager, A. and Zolotarev, E.},
        title = {Counter-examples to a conjecture of Karpenko via truncated Brown-Peterson cohomology},
        year = {2026}
    }
