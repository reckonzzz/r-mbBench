# mbBench for R

Load pre-generated microbiome benchmarks, evaluate method results and generate new profiles. Requires R 4.1 or later.

## Install

```r
install.packages("remotes")  # Run once if needed
remotes::install_github("reckonzzz/r-mbBench")
library(mbBench)
```

## Quick start

**Evaluate a method.** Load a [pre-generated benchmark](https://doi.org/10.5281/zenodo.24076644) directly in R. Run your method on `data$counts`, then score its results with `mb_evaluate_da()` or `mb_evaluate_batch()`.

```r
data <- mb_dataset(preset = "da-null-balanced", seed = 42L)
```

**Generate new profiles.** Configure the Python backend once before using a pretrained model or logit bank.

```r
install.packages("reticulate")  # Run once if needed
mb_install_backend()
```

To use a local [Python checkout](https://github.com/reckonzzz/mbBench), call `mb_install_backend(source = "../mbBench")` from this repository with the two checkouts side by side.

Choose `ae` (default), `vae`, or `ae-system` from the [pretrained model release](https://doi.org/10.5281/zenodo.24080099). The `ae-system` model adds Body-System conditioning.

```r
model <- mb_download_model("ae")
result <- mb_generate(model = model, n_samples = 20L, seed = 42L)
```

**Reuse stored logits.** Use the [core logit bank](https://doi.org/10.5281/zenodo.24082591) to generate profiles on CPU without diffusion sampling.

```r
bank <- mb_download_bank()
result <- mb_generate(bank = bank, n_samples = 20L, seed = 42L)
```

Use `mb_benchmark()` to create a benchmark and `help(package = "mbBench")` for function documentation.

## R functions

Use `?mb_dataset` for function help or `help(package = "mbBench")` for the package index.

| R function | Purpose |
| --- | --- |
| `mb_install_backend()` | Install or configure the Python backend used for simulation. |
| `mb_backend_status()` | Report the Python backend configuration and status. |
| `mb_list_datasets()` | List the available pre-generated benchmark scenarios. |
| `mb_show_dataset()` | Show the available seeds and catalog details for a scenario. |
| `mb_download_dataset()` | Download and load an official pre-generated benchmark dataset. |
| `mb_dataset()` | Load a benchmark by scenario and seed, from a ZIP, or from a directory. |
| `mb_download_model()` | Download and load an official pretrained model. |
| `mb_load_model()` | Load an installed official model or a custom model directory. |
| `mb_download_bank()` | Download and load a pre-generated logit bank. |
| `mb_load_bank()` | Load a logit bank from an extracted directory. |
| `mb_list_conditions()` | List accepted condition values for a loaded model or bank. |
| `mb_generate()` | Generate relative abundances or observed counts from a model or bank. |
| `mb_generate_from_condition_table()` | Generate samples for each row of a condition table. |
| `mb_list_presets()` | List built-in benchmark scenarios. |
| `mb_show_preset()` | Inspect the complete configuration for a preset. |
| `mb_list_parameters()` | Show benchmark parameters, defaults and accepted values. |
| `mb_write_config_template()` | Write a YAML template for a custom benchmark. |
| `mb_benchmark()` | Generate a benchmark with known truth from a preset or custom configuration. |
| `mb_da_result()` | Build a complete differential-abundance decision table from significant and tested taxa. |
| `mb_evaluate_da()` | Score differential-abundance decisions against benchmark truth. |
| `mb_evaluate_da_many()` | Evaluate a manifest of differential-abundance results and summarize repetitions. |
| `mb_summarize_da()` | Aggregate differential-abundance metrics across repetitions. |
| `mb_evaluate_batch()` | Score corrected abundance matrices against a benchmark. |
| `mb_summarize_batch()` | Aggregate batch-correction metrics and scores across repetitions. |
