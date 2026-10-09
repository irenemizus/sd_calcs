# sd_calcs

A script for finding related pairs of calculated and observed rovibrational
energy levels of a molecule and obtaining their deviations and standard
deviation (sd) values.

For each experimental level the program locates the best-matching calculated
level **within the same `(J, symmetry)` group** by minimising the
observed − calculated energy difference. Every matched pair is then labelled
according to an energy threshold `eps`:

- **`FOUND`** — `|obs − calc| <= eps`
- **`OUTLIER`** — `|obs − calc| > eps`

From these labelled pairs the program computes the standard deviation of the
deviations, separately for low-, mid- and high-energy ranges, and can also
produce fit-input files.

## Supported molecules

- `H2-16O` (water)
- `N2O` (nitrous oxide)
- `N2O-556` (an `N2O` isotopologue; RITZ experimental data only)

## Supported data

- **Molecules:** `H2-16O`, `N2O`, `N2O-556`
- **Files with calculated data:** `fort.14`-like, ExoMol states-like, HITRAN
- **Files with observed data:** MARVEL-like, RITZ-like (for `N2O-556`)
- **Output formats:** Yurchenko's fit-input file, comparison files, sd summary files
- **Output labelling:**
  - symmetries `A1, A2, B1, B2` (plus `UKN` for undetermined)
  - 6 quantum numbers for `H2-16O` (`v1, v2, v3, J, Ka, Kc`)
  - 4 quantum numbers for `N2O` (`v1, v2, J, l`)
  - 5 quantum numbers for `N2O-556` (`v1, v2, v3, J, l`)

## Requirements

- Python 3 (tested with 3.8)
- Standard library only — no third-party packages required.

## Repository layout

```
.
├── main.py               # main entry point: obs-vs-calc matching + sd
├── obs-calc_comp.py      # add-on: compare previous (HITRAN/ExoMol) calc data vs the current comparison
├── states.py             # data model (State, States, quantum numbers, symmetries, statuses, comparison lists)
├── formats.py            # input-file parsers (fort.14, ExoMol, HITRAN, MARVEL exp, RITZ exp)
├── README.md             # this file
├── input/                # input data, one sub-directory per molecule (git-ignored)
│   ├── H2-16O/           #   fort.14-* files + observed-data file
│   ├── N2O/
│   └── N2O-556/          #   fort.14-* files + RITZ observed-data file
└── output/               # generated results, one sub-directory per molecule (git-ignored)
```

> `input/` and `output/` are not tracked in git. The program reads from
> `input/<mol_name>/` and writes to `output/<mol_name>/`.

## How the matching works

For a given `J` the calculated levels are grouped by symmetry. For each
experimental state the program:

1. takes the first calculated level with the same `J` and symmetry,
2. walks forward through the remaining same-group levels, keeping the one with
   the smallest `|obs − calc|`, and
3. stops once the sign of the discrepancy changes or the absolute deviation
   starts increasing.

The winning calculated level is then removed so it cannot be re-matched, and the
experimental state is tagged `FOUND` or `OUTLIER` using the `eps` threshold.

### Standard-deviation ranges

The sd of the deviations is computed over three energy bands, split by the
experimental energy `E_exp`:

- **low** — `E_exp <= E_tr_l`
- **middle** — `E_tr_l < E_exp < E_tr_h`
- **high** — `E_exp >= E_tr_h`

For each band both the sd over `FOUND` states and the sd over `OUTLIER` states
are reported.

## Usage

### `main.py`

```
python main.py [options]
```

| Option | Description | Default |
|--------|-------------|---------|
| `--mode` | `all` (process every `fort.14*` file in the input folder) or the name of a single file (also handles a single ExoMol states file) | `all` |
| `--mol_name` | molecule name: `H2-16O`, `N2O`, or `N2O-556` | `H2-16O` |
| `--file_exp_name` | name of the observed-data file inside `input/<mol_name>/` | — |
| `--out_file_exp_name` | base name for the Yurchenko fit-input output file (extended automatically) | `ens_Yur_format.txt` |
| `--out_file_comp_name` | base name for the full comparison output file (extended automatically) | `comp.txt` |
| `--E_zero` | zero of energy for the calculated data | — |
| `--E_tr_l` | upper edge of the low-energy band (cm⁻¹) | `15000.0` |
| `--E_tr_h` | upper edge of the middle-energy band (cm⁻¹) | `25000.0` |
| `--eps` | deviation threshold (cm⁻¹) separating `FOUND` from `OUTLIER` | `0.2` |
| `--Jmax` | maximum `J` when using ExoMol states-like calculated data | `100` |
| `--make_comp_files` | `True` (generate comparison/fit-input files) or `False` (reuse existing ones from the output folder) | `True` |
| `--N_out_rel` | percentage of the worst-predicted levels to treat as outliers for an extra all-`J` sd calculation | `0.0` |

#### Example

```bash
# Compare all fort.14* files for water against the MARVEL observed levels,
# eps = 0.2 cm-1, and generate the comparison + fit-input files.
python main.py \
    --mode all \
    --mol_name H2-16O \
    --file_exp_name W2024_energy_levels_h2-16o.txt \
    --eps 0.2
```

```bash
# Compare a single ExoMol states file (all J up to --Jmax).
python main.py \
    --mode H2O_J69_all.states \
    --mol_name H2-16O \
    --file_exp_name W2024_energy_levels_h2-16o.txt \
    --Jmax 69 \
    --eps 0.2
```

```bash
# N2O against observed levels, and also mark the worst 10 % as outliers.
python main.py \
    --mode all \
    --mol_name N2O \
    --file_exp_name <n2o_obs_file> \
    --N_out_rel 10.0
```

```bash
# N2O-556 isotopologue against RITZ observed levels (J is the first column in the input file).
python main.py \
    --mode all \
    --mol_name N2O-556 \
    --file_exp_name RITZ.556.levels.v2 \
    --eps 0.2
```

### Outputs of `main.py`

Written to `output/<mol_name>/`:

- **`<out_file_comp_name>...`** — the full comparison (one line per matched
  state: `J`, symmetry, `N`, `E_exp`, `E_calc`, `E_diff`, quantum numbers,
  weight, status).
- **`...+sd...`** — the same comparison with the per-band sd values appended.
- **`<out_file_exp_name>...`** — Yurchenko fit-input file.
- **`out_file_all_sd_*.txt`** — sd values for every `J` and for all `J`s
  combined.
- When `--N_out_rel > 0`: extra files with the worst levels re-marked as
  outliers and the corresponding sd values.
- When `--make_comp_files False`: a `..._valid` file holding only the
  zero-weighted states and their sd.

> Changing `--eps` requires re-running with `--make_comp_files True`; every new
> `eps` value needs a recalculation from the beginning.

### `obs-calc_comp.py`

A small add-on to `main.py`. Given an already-generated "current calculated
vs observed" comparison file and a *previous* calculated dataset
(HITRAN / ExoMol states), it re-matches the previous calculated levels to each
compared state (by quantum numbers, or by best fit) and reports how the
deviations changed relative to the previous calculation.

```
python obs-calc_comp.py [options]
```

| Option | Description | Default |
|--------|-------------|---------|
| `--mol_name` | molecule name: `H2-16O` or `N2O` | `H2-16O` |
| `--file_comp_name` | name of the pre-generated "current calc vs obs" comparison file | — |
| `--file_old_calc_name` | name of the previous (HITRAN/ExoMol) calculated-data file | — |
| `--out_file_comp_name` | base name for the output comparison file | `comp.txt` |
| `--eps` | deviation threshold used in the best-fit matching (cm⁻¹) | `0.2` |
| `--Jmax` | maximum `J` for comparison | `100` |

It reads from and writes to
`input/<mol_name>_old_calc_comp/` and
`output/<mol_name>_old_calc_comp/` respectively.

Each compared state is annotated with a mark indicating the change of the
deviation relative to the previous calculation: `≈` for a relative change
between 0.5 and 1.0, `!!!` for a relative change greater than 1.0.
