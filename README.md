# Perplexity

## Usage

```
USAGE:
    perplexity eval [FLAGS] [OPTIONS] --output <output> --quants <train_ecs> --test_ecs <test_salmon_dir> --train_ecs <train_ecs>

FLAGS:
    -h, --help             Prints help information
    -t, --smooth_per_tx    use reads per tx smoothing
    -V, --version          Prints version information

OPTIONS:
    -o, --output <output>               output results path
        --quants <train_ecs>            Training set quants
        --smoothing <smoothing>         TPM based smoothing [default: 1e-8]
        --test_ecs <test_salmon_dir>    Test set equivalence class counts
        --train_ecs <train_ecs>         Training set equivalence class counts
```

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
