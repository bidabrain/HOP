# HOP: A Group-Finding Algorithm for N-Body Simulations

HOP is an algorithm for identifying particle groups (halos) in N-body simulations based on local density estimates. It is spatially adaptive, coordinate-free, and numerically straightforward.

This repository is a fork of the original HOP distribution from the [CfA HOP page](https://lweb.cfa.harvard.edu/~deisenst/hop/).

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
