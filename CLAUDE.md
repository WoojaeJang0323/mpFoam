# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`mpFoam` is a custom OpenFOAM solver for mineral precipitation modeling. It extends `icoReactingMultiphaseInterFoam` using a VOF (Volume of Fluid) approach with two phases (liquid and solid) and a custom precipitation mass transfer model. The active solver is in `solver/MVP1/`. OpenFoam v1912 (ESI) is currently used for this solver.

## Build Commands

Requires a sourced OpenFOAM environment (`source /path/to/OpenFOAM/etc/bashrc`).

```sh
# Build everything (from solver/MVP1/)
cd solver/MVP1
./Allwmake

# Clean build artifacts
./Allwclean
```

`Allwmake` calls `wmake` in this order:
1. `phasesSystem` → library `IncompressibleMultiphaseSystemsNew`
2. `massTransferModels` → library `massTransferModelsNew`
3. `CompressibleMultiPhaseTurbulenceModels`
4. Solver executable `mpFoam`

Libraries install to `$FOAM_USER_LIBBIN`; the executable to `$FOAM_USER_APPBIN`.

## Running Cases

Each case directory has two run scripts:

```sh
# Parallel run (decomposes, runs MPI, reconstructs)
./Allrun

# Single-process run (no decomposition)
./singRun
```

The general workflow: `blockMesh` (or copy pre-built mesh) → restore `.orig` initial conditions → `setFields` → `decomposePar` → `mpirun mpFoam` → `reconstructPar`.

Cases with complex geometry (`precFlowCylinder`, `precFlowSphere`) use `snappyHexMesh` for meshing (a separate `case/meshing/` step produces the mesh, which is then copied into `constant/polyMesh/`).

## Solver Architecture

### Main time-step loop (`solver/MVP1/mpFoam.C`)
1. `fluid.solve()` — solves phase fraction equations (CMULES algorithm, VOF)
2. PIMPLE loop:
   - `UEqn.H` — momentum equation
   - `CEqn.H` — concentration transport equation
   - `fluid.correct()` — interface compression and mass transfer source terms
   - `pEqn.H` — pressure correction

### Custom Libraries

**`phasesSystem/`** — Multiphase system framework:
- `phaseSystem`: Base class managing all phases and their interactions
- `multiphaseSystem`: Adds CMULES VOF solver, mass transfer source coupling, mass balance tracking, CSV output
- `MassTransferPhaseSystem`: Connects mass transfer models to the phase system
- `phaseModel/`: Phase types — `MovingPhaseModel` (liquid), `StaticPhaseModel` (solid)
- `interfaceModels/porousModels/VollerPrakash`: Darcy-like drag term to suppress flow inside solid (parameters `Cu`, `solidPhase` set in `phaseProperties`)

**`massTransferModels/precipitate/`** — Core precipitation physics:
- `Kexp()` method computes the mass flux from liquid → solid
- Surface area density is computed as `2 * |∇alpha.liquid|` (gradient method) or a uniform constant (`smoothSurface` option)
- A `Cmask` field marks cells adjacent to existing solid, restricting where precipitation occurs
- Rate is capped each timestep to prevent the concentration from going negative

### Concentration Equation (`CEqn.H`)
Diffusivity is active only in the liquid phase: `diffEff = diffInLiquid × alpha_liquid²`. This concentrates diffusion where there is liquid and naturally zeroes it out in solid regions.

## Key Configuration Files

**`constant/inputParameters`** (per-case):
```
molW          18.0;      // Molar weight of solid [g/mol]
diffInLiquid  1e-6;      // Diffusion coefficient in liquid [m²/s]
```

**`constant/phaseProperties`** (per-case):
```
massTransferModel
(
    (liquid to solid)
    {
        type            precipitate;
        C               1E-3;       // Reaction rate [m/s]
        Tactivate       310;        // Activation temperature [K]
        Cactivate       1;          // Activation concentration [mol/m³]
        Mv              18;         // Molar weight [g/mol]
        alphaSolidMin   0.9;        // Threshold to detect solid cells
        smoothSurface   false;      // true = constant area density, false = gradient-based
        smoothAreaDensity 1.0;      // Used only when smoothSurface is true [1/m]
    }
);
```

## Fields

| Field | Description |
|-------|-------------|
| `alpha.liquid` / `alpha.solid` | Phase volume fractions (sum ≤ 1) |
| `C` | Solute concentration [mol/m³] |
| `diffCo` | Effective diffusion coefficient (computed from `alpha.liquid` each step) |
| `areaDensityGrad` | Interface area density [1/m], written at output times |
| `Cmask` | Binary field marking cells eligible for precipitation |
| `SuOut` | Source term output field from the precipitation model |
| `dmdt.liquidToSolid` | Mass transfer rate field |

## Solver Versions

| Directory | Status |
|-----------|--------|
| `solver/v0.1`, `v0.2`, `v0.3` | Earlier development versions, not actively used |
| `solver/MVP1` | Current version; includes precipitation model and validation cases |
