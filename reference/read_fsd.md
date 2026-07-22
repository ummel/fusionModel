# Read Fusion Output Files From Disk

`read_fsd()` efficiently loads microdata fusion outputs (`.fsd` files)
or standard Fast Storage (`.fst`) datasets into memory. Because fusion
files can become extremely large when generating multiple implicates,
`read_fsd()` leverages partial reading capabilities to load only the
specific rows, columns, or implicates needed, drastically reducing
memory footprint and processing time.

## Usage

``` r
read_fsd(
  path,
  columns = NULL,
  M = 1,
  df = NULL,
  cores = max(1, parallel::detectCores(logical = FALSE) - 1)
)
```

## Arguments

- path:

  Character. Path to a valid `.fsd` or `.fst` file on disk, typically
  created by calling
  [`fuse`](https://ummel.github.io/fusionModel/reference/fuse.md) with
  the `fsd` argument specified.

- columns:

  Character vector. Specific column names to read from disk. The default
  (`NULL`) loads all available columns in the dataset.

- M:

  Numeric or Integer. The maximum number of implicates to read (e.g.,
  `M = 5` loads implicates `1` through `5`). Set `M = Inf` to load all
  implicates available in the file. Default is `1`. Ignored if the
  target dataset does not contain an `M` column.

- df:

  Data frame or data.table. Optional filtering subset. If provided, only
  rows in the file matching the combination of key values in `df` will
  be returned. Default is `NULL`.

- cores:

  Integer. Number of CPU threads to assign to `fst` for multi-threaded
  decompression. Defaults to available logical cores minus one (minimum
  1).

## Value

A [`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html)
containing the requested columns and rows. Original `data.table` key
attributes are automatically preserved if present in the source file.

## Details

The `.fsd` format (Fusion Data file) is functionally identical to the
high-performance binary format supplied by the fst package, optimized
for fast random access and compression.

**Subsetting Mechanics:**

- **Implicate Filtering:** When the dataset contains multiple implicates
  (identified by column `M`), `read_fsd()` calculates contiguous row
  indices corresponding to `M <= max_M` before reading, preventing
  unnecessary disk I/O for unused implicates.

- **Row Subsetting (`df`):** If a matching dataset `df` is provided:

  - For files **\< 100 MB**, the relevant implicate range is fully
    loaded into memory and filtered using an efficient inner join via
    [`join`](https://fastverse.org/collapse/reference/join.html).

  - For files **\>= 100 MB**, low-memory index matching is performed
    out-of-core using high-speed C-level matching via
    [`fmatch`](https://fastverse.org/collapse/reference/fmatch.html)
    (`%iin%`), loading only matching record indices.

## Examples

``` r
if (FALSE) { # \dontrun{
# Build a fusion model using RECS microdata
fusion.vars <- c("electricity", "natural_gas", "aircon")
predictor.vars <- names(recs)[2:12]
fsn.path <- train(data = recs, y = fusion.vars, x = predictor.vars)

# Write fusion implicates directly to disk during fusion
recipient <- recs[predictor.vars]
fsd_file <- file.path(tempdir(), "results.fsd")
fuse(data = recipient, fsn = fsn.path, M = 3, fsd = fsd_file)

# Read only the first implicate (M = 1) for all variables
sim_m1 <- read_fsd(path = fsd_file, M = 1)
head(sim_m1)

# Read specific columns across all 3 implicates
sim_sub <- read_fsd(
  path = fsd_file,
  columns = c("M", "electricity"),
  M = 3
)
head(sim_sub)
} # }
```
