# Two-Inlet Mixing Precipitation

This run case injects anion and cation solutions from two separate inlet patches:

- `inletAnion`: `Canion = 8 mol/m^3`, `Ccation = 0`
- `inletCation`: `Ccation = 8 mol/m^3`, `Canion = 0`

The two streams move through a 2D channel and mix along the midline. A small
solid seed is initialized at `x = 2..3 mm`, `z = 3.6..4.4 mm`
(`0.002..0.003 m`, `0.0036..0.0044 m` in `setFieldsDict`) so
`mixingPrecipitation` has an initial nucleation surface inside the mixing zone.
The top and bottom wall patch `walls` is also listed in `reactiveWallPatches`.

Run with:

```bash
source /usr/lib/openfoam/openfoam1912/etc/bashrc
cd /home/wjj/mpFoam/runs/twoInletMixingPrecipitation
./Allrun
```

Important settings:

- `constant/inputParameters`: enables `useTwoSpecies yes`
- `constant/phaseProperties`: uses `type mixingPrecipitation`
- `system/blockMeshDict`: split left boundary into `inletAnion` and `inletCation`
