# HOP: A Group-Finding Algorithm for N-Body Simulations

HOP is an algorithm for identifying particle groups (halos) in N-body simulations based on local density estimates. It is spatially adaptive, coordinate-free, and numerically straightforward.

This repository is a fork of the original HOP distribution from the [CfA HOP page](https://lweb.cfa.harvard.edu/~deisenst/hop/). Versions 1.0 and 1.1 are the work of the original authors. **All modifications from v1.2 onward are by [@bidabrain](https://github.com/bidabrain)** and are not associated with the original authors.

## Authors

**Daniel J. Eisenstein** and **Piet Hut**

- Daniel J. Eisenstein — deisenstein@cfa.harvard.edu
- Piet Hut — piet@sns.ias.edu

The algorithm is described in full in:

> Eisenstein, D. J. & Hut, P. (1998). "HOP: A New Group-Finding Algorithm for N-Body Simulations." *Astrophysical Journal*, 498, 137. (Published May 1, 1998)

The kd-tree search engine (`smooth.c`, `kd.c`) was written by **Joachim Stadel** and the **NASA HPCC ESS at the University of Washington Department of Astronomy**, distributed as part of SMOOTH v2.01.

## How HOP Works

HOP assigns a local density to every particle via an adaptive smoothing kernel, then traces density-increasing paths to group particles together:

1. **Density estimation** — Compute the density around each particle using an adaptive kernel scaled to the distance of the N_dens nearest neighbor.
2. **Hopping** — Each particle records the ID of its densest neighbor among N_hop neighbors. Starting from every particle, hop along these links until reaching a local density maximum; all particles converging to the same maximum form one group.
3. **Boundary cataloging** — For each pair of adjacent groups, find the highest-density boundary particle pair and record the average density of that boundary.
4. **Group merging** — Viable groups (peak density > delta_peak) are merged across shared boundaries with density > delta_saddle. Unviable groups are absorbed into their best viable neighbor.
5. **Outer density cut** — Particles with density below delta_outer are excluded from all groups.

Steps 1–4 are performed by `hop`, and steps 5–6 by `regroup`, allowing the density thresholds to be varied cheaply without re-running the expensive neighbor searches.

## Source Files

| File | Description |
|---|---|
| `hop.c` | Main program: density estimation, hopping, boundary cataloging |
| `hop_input.c` | **User-customizable** input routine for reading simulation files |
| `smooth.c` | kd-tree search engine (from SMOOTH v2.01, Stadel/HPCC) |
| `smooth.h` | Header for smooth.c |
| `kd.c` | kd-tree construction utilities (from SMOOTH v2.01) |
| `kd.h` | Header for kd.c |
| `regroup.c` | Group merging and density threshold application |
| `slice.c` | Utility routines for managing simulation data |
| `slice.h` | Header for slice.c |
| `Makefile` | Build file |

## Compilation

1. **Customize the input routine.** Edit `hop_input.c` to read your simulation file format. Decide whether all particles have equal mass or not.

2. **Edit the Makefile.** Set `CC` to your compiler, `CFLAGS` to your optimization flags, and `LIBS` if you need to replace `libm`. If particles have unequal masses, uncomment the `-DDIFFERENT_MASSES` line.

3. **Build.**

   ```sh
   make
   ```

   This produces two executables: `hop` and `regroup`. Run `make clean` to remove intermediate object files.

## Usage

### Step 1 — Run `hop`

```sh
hop -in <simulation_file> -p <box_size> -nd <N_dens> -nh <N_hop> -nm <N_merge> -dt <density_threshold> -out <output_root>
```

Key flags:

| Flag | Default | Description |
|---|---|---|
| `-in file` | stdin | Input simulation file |
| `-p float` | none | Periodic box size (required for cosmological sims) |
| `-nd int` | 64 | N_dens: neighbors used for density smoothing |
| `-nh int` | 16 | N_hop: neighbors searched when hopping |
| `-nm int` | 4 | N_merge: neighbors searched for group boundaries (must be < N_hop) |
| `-dt float` | none | Density threshold; particles below this are excluded from `.gbound` |
| `-den file` | — | Read pre-computed densities from file (skips density step) |
| `-out file` | — | Root name for output files |
| `-densityonly` | — | Compute densities only, then quit |

Output files: `<root>.den` (densities), `<root>.hop` (raw group memberships), `<root>.gbound` (group boundary info).

**Example:**

```sh
hop -in my_sim -p 1 -nd 48 -nh 20 -nm 5 -dt 1 -out my_hop_out
```

To re-run with a different N_hop while reusing saved densities:

```sh
hop -in my_sim -p 1 -den my_hop_out.den -nh 16 -nm 5 -dt 1 -out my_hop_out2
```

### Step 2 — Run `regroup`

```sh
regroup -root <hop_output_root> -douter <delta_outer> -dsaddle <delta_saddle> -dpeak <delta_peak> -mingroup <N> -out <output_root>
```

Key flags:

| Flag | Default | Description |
|---|---|---|
| `-root file` | — | Common root for `.hop`, `.den`, `.gbound` inputs |
| `-douter float` | — | delta_outer: minimum density for group membership |
| `-dsaddle float` | — | delta_saddle: minimum boundary density to merge two viable groups |
| `-dpeak float` | — | delta_peak: minimum peak density for a group to be viable |
| `-mingroup int` | 10 | Minimum group size to output |
| `-pipe` | — | Write `.tag` to stdout instead of disk |
| `-f77` | — | Write `.tag` in FORTRAN unformatted binary format |

Output files: `<root>.tag` (final group memberships), `<root>.size` (group sizes), `<root>.gmerge` (merge log).

**Example:**

```sh
regroup -root my_hop_out -douter 80 -dsaddle 140 -dpeak 160 -mingroup 8 -out my_final_catalog
```

## Recommended Parameters

From the authors' paper, these parameters gave robust results on cosmological simulations:

- N_dens = 64, N_hop = 16, N_merge = 4
- delta_outer = delta_saddle / 2.5 = delta_peak / 3

## Version History

- **v1.0** — Initial release
- **v1.1** (April 2, 2002) — Fixed a bug in `regroup.c` that could cause a crash (not data corruption) for very large particle numbers
- **v1.2** (June 7, 2026) — Modifications by [@bidabrain](https://github.com/bidabrain). All changes from this version onward are independent of the original authors. Bug fixes in `regroup.c` and `hop.c`:
  - **`regroup.c` — crash**: `sort_groups()` called `fclose(f)` unconditionally; when no `.size` output file was requested, `f` was never initialized, causing undefined behavior / crash
  - **`regroup.c` — stdout corruption**: `writetags()` and `writetagsf77()` called `fclose(f)` even when writing to stdout (pipe mode), closing stdout and corrupting all subsequent output
  - **`regroup.c` — silent data truncation**: `readtags()` did not check the return value of `fread` for the tag array; a truncated `.hop` file would silently assign wrong group IDs to all particles after the cut point
  - **`regroup.c` — uninitialized density data**: `merge_groups_boundaries()` did not verify that all `ngroups` entries were read before the `###` separator; a premature `###` left `gdensity[]` entries uninitialized, causing wrong peak-density decisions
  - **`regroup.c` — missing header read check**: `densitycut()` did not check the return value of `fread` for the particle count header
  - **`regroup.c` — buffer overflow**: filename buffers in `parsecommandline()` were `malloc(80)`; root paths longer than ~75 characters would overflow; increased to `malloc(1024)`
  - **`regroup.c` — line buffer too small**: `merge_groups_boundaries()` used `char line[80]`; `.gbound` group lines with large integer IDs or float values could exceed 79 characters; increased to `char line[256]`
  - **`regroup.c` / `hop.c` — binary files in text mode**: all binary input/output files (`.hop`, `.den`) were opened with `"r"`/`"w"` instead of `"rb"`/`"wb"`; results are correct on Unix but wrong on Windows
  - **`regroup.c` / `hop.c`**: `void main()` changed to `int main()` (non-standard return type)
  - **`regroup.c`**: added `#include <float.h>` for `FLT_MAX`
