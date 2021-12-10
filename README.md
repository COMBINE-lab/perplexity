# Perplexity

## Usage

```
USAGE:
    perplexity eval [OPTIONS] --output <output> --quants <train_ecs> --test_ecs <test_salmon_dir> --train_ecs <train_ecs>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -o, --output <output>               Write results to output
        --quants <train_ecs>            Training set quant.sf
        --smoothing <smoothing>         Smoothing parameter for various strategies [default: 1e-8]
    -m, --smoothing_strategy <mode>     One of {TX, TPM, LGT} [default: LGT]
        --test_ecs <test_salmon_dir>    Test set equivalence class file
        --train_ecs <train_ecs>         Training set equivalence class file
```

Use `perplexity --help` for more options

## Snakemake workflow Usage

See `snakefiles/config.yml` for details regarding input filepaths and conventions for setting up k-fold quantify-then-validate procedure. Supplied `perplexity.snk` snakefile runs Salmon on k-folds then computes perplexity of held-out reads on quantified reads. `kfold.snk` is used to set up quantify-validate fragment sets.

### Perplexity of Salmon abundance estimates across K-folds
```
snakemake  --snakefile snakefiles/perplexity.snk --configfile <config file> -j <n jobs>
```
  
### Set up train-test split
```
snakemake --snakefile snakefiles/kfold.snk --configfile <config file> -j <n jobs>
```

## Build

```
cargo build --release
```
