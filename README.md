# glmGamPoi\[2\] <img src="logo.svg" alt="logo" height="128" align="right"/>

> Fit Gamma-Poisson Generalized Linear Models Reliably.

> [!IMPORTANT]
> this repo is a fork of upstream [glmGamPoi](https://github.com/const-ae/glmGamPoi).
> it is specifically designed to optimize memory use and runtime, without sacrificing the library's existing accuracy or features.
> more details regarding using the fork, its scope, and relevant modifications are discussed in [fork.md](./fork.md).

Pronounciation: [`dʒi əl əm ɡam ˈpwɑ`](http://ipa-reader.com/?text=d%CA%92i%20%C9%99l%20%C9%99m%20%C9%A1am%20%CB%88pw%C9%91)

The core design aims of `glmGamPoi[2]` are:

- Fit Gamma-Poisson models on arbitrarily large or small datasets
- Be faster and more memory efficient than alternative methods, such as `DESeq2` or `edgeR`
- Calculate exact or approximate results based on user preference
- Support in-memory or on-disk data
- Follow established conventions around tools for RNA-seq analysis
- Present a simple user-interface
- Avoid unnecessary dependencies
- Make integration into other tools easy

## Installation

You can build and install [`glmGamPoi2`](https://github.com/goeva-lab/glmGamPoi) from the automatically automatically built package archive published from the latest [GitHub release](https://github.com/goeva-lab/glmGamPoi/releases/latest/download/glmGamPoi.tar.gz).

This can be done using the [`pak`](https://pak.r-lib.org/) (for installation) R package via the following snippet:

```R
pak::pak("url::https://github.com/goeva-lab/glmGamPoi/releases/latest/download/glmGamPoi2.tar.gz")
```

<details><summary>Details regarding optimized build instructions</summary>

To build and install the package with the maximal level of optimizations enabled for your particular machine, the  [`withr`](https://withr.r-lib.org/) package can be used as follows:

```R
shared_flags <- "-fuse-ld=lld -flto -ffat-lto-objects -fopenmp"
withr::with_makevars(assignment = "=", c(
    CXX17STD = "-std=c++17",
    CXX17FLAGS = sprintf("%s -O3 -march=native", shared_flags),
    LDFLAGS = sprintf("%s -L/usr/lib64/R/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,-z,pack-relative-relocs", shared_flags)
), {
    pak::pak("url::https://github.com/goeva-lab/glmGamPoi/releases/latest/download/glmGamPoi2.tar.gz")
})
```

</details>

If you use this package in a scientific publication, please make sure to cite the [original glmGamPoi publication](https://doi.org/10.1093/bioinformatics/btaa1009):

> glmGamPoi: Fitting Gamma-Poisson Generalized Linear Models on Single
> Cell Count Data  
> Constantin Ahlmann-Eltze, Wolfgang Huber  
> Bioinformatics; 2020-12-09; doi:
> <https://doi.org/10.1093/bioinformatics/btaa1009>

## Example

An example of library use can be seen below:

```R
library(glmGamPoi2)

counts <- rnbinom(n = 10, mu = 5, size = 1/0.7)

# design = ~ 1 means that an intercept-only model is fit
fit <- glm_gp(counts, design = ~ 1)
```

## Acknowlegments

Original work for the [original glmGamPoi package](https://github.com/const-ae/glmGamPoi) was supported by the EMBL International PhD Programme and the European Research Council Synergy grant DECODE under grant agreement No. 810296.

Work for the [glmGamPoi2 fork](https://github.com/goeva-lab/glmGamPoi2) was supported by startup funds from the University of Toronto under the supervision of Dr. Aleksandrina Goeva.
