# Figure 4 Validation Case

This case reproduces the `mpFoam` side of Figure 4 from:

- Yang, Stack, Starchenko (2021), *Micro-continuum approach for mineral precipitation*.

Implemented setup:

- domain size: `2.1 mm × 1.5 mm × 1.5 mm`
- sphere radius: `100 um`
- sphere center: `x = 0.7 mm`, `y = 0`, `z = 0`
- inlet velocity: `1e-4 m/s` (`100 um/s`)
- equal cation/anion concentrations: `0.32 mol/m^3`
- saturation concentration: `Cs = 0.0105 mol/m^3`
- diffusion coefficient: `1.4e-9 m^2/s`
- kinematic viscosity: `0.89e-6 m^2/s`

Model choices in this repo:

- uses `mixingPrecipitation`
- uses two equal species fields `Ccation` and `Canion` to represent the paper's
  equal-ion assumption in each cell
- uses concentration-only saturation:
  `Omega = Ccation * Canion / Ksp`
- `Ksp = Cs^2 = 1.1025e-4`
- `nPrec = 1.0`

Important implementation note:

- The current wall-reactive model in this repo does not yet compute wall area
  density from boundary face area directly, so this case uses
  `smoothSurface true` with `smoothAreaDensity 1e5` to activate precipitation on
  the `sphere` wall patch.

Run with:

```bash
source /usr/lib/openfoam/openfoam1912/etc/bashrc
cd /home/wjj/mpFoam/runs/figure4ValidationEq7
./Allrun
```

There is no ALE solver implementation in this workspace, so only the `mpFoam`
validation-side case is provided here.
