# LMX_data

> **A Database for Design and Optimization of Halide Based Solid Electrolytes for LIBs**

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![ASE](https://img.shields.io/badge/ASE-3.x-green.svg)](https://wiki.fysik.dtu.dk/ase/)
[![CLEASE](https://img.shields.io/badge/CLEASE-latest-orange.svg)](https://clease.readthedocs.io/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Table of Contents

- [Overview](#overview)
- [Background](#background)
- [Repository Structure](#repository-structure)
- [Crystal Structure Templates](#crystal-structure-templates)
- [Script Descriptions](#script-descriptions)
- [Installation & Dependencies](#installation--dependencies)
- [Usage](#usage)
- [CrysGen Class](#crysgen-class-api-reference)
- [Template Reference Table](#template-reference-table)
- [Notes & Known Issues](#notes--known-issues)
- [Author](#author)

---

## Overview

Generating Li-M-X (Lithium–Metal–Halide) solid electrolytes datset using **Cluster Expansion (CE)** methodology via the [CLEASE](https://clease.readthedocs.io/) library and [ASE](https://wiki.fysik.dtu.dk/ase/) (Atomic Simulation Environment).

The toolkit supports **six distinct crystal symmetry prototypes** (space groups 11, 12, 62, 164, and 165) across four generation modes:

| Mode | Purpose |
|---|---|
| `GenRandom` | Mixed metal + halogen substitution at fixed halogen concentration |
| `GenRandomAll` | Multi-generation mixed metal + halogen substitution |
| `GenRandomMetal` | Metal-site substitution only, at a fixed metal concentration |
| `GenRandomDefect` | Defect structures with Li vacancies and anti-site defects |
| `GenRandomDefectAll` | Multi-generation defect structure generation (all templates) |

---

## Background

Lithium halide solid electrolytes of the general formula **Li-M-X** (where M = trivalent metal, X = halogen) are promising candidates for all-solid-state lithium-ion batteries due to their high ionic conductivity and electrochemical stability. Systematic exploration of composition and defect spaces is critical to identifying optimal materials.
---

## Repository Structure

LMX_data/

|--- GenRandom.py           # Single-generation mixed metal+halogen CE structures

|--- GenRandomAll.py        # Multi-generation mixed metal+halogen CE structures

|--- GenRandomMetal.py      # Single-generation metal-only CE structures

|--- GenRandomDefect.py     # Single-generation defect (vacancy/anti-site) CE structures

|--- GenRandomDefectAll.py  # Multi-generation defect CE structures (all templates)

|--- README.md              # This file

> **Note:** CIF template files are expected at a local path (e.g., `LMX_prototypes`). You must update these paths to match your local setup before running any script.

---

## Crystal Structure Templates

Six prototype crystal structures (named after their space group and halide type) are supported:

| Template Name | Space Group | Space Group Symbol | Prototype Compound | Primary Halide |
|---|---|---|---|---|
| `F11` | 11 | P 1 2₁/m 1 | LiInF₄ | Fluoride |
| `Cl12` | 12 | C 1 2/m 1 | Li₃ScCl₆ | Chloride |
| `Br12` | 12 | C 1 2/m 1 | Li₃ErBr₆ | Bromide |
| `F15` | 15 | C 1 2/c 1 | Li₃GaF₆ | Fluoride |
| `Cl62` | 62 | P n m a | Li₃YbCl₆ | Chloride |
| `Cl164` | 164 | P -3 m 1 | Li₃ErCl₆ | Chloride |
| `F165` | 165 | P -3 c 1 | Li₃ScF₆ | Fluoride |

Each template defines:
- **Lattice parameters** (`a, b, c, α, β, γ`) read from a reference CIF file
- **Wyckoff positions** for all sublattice sites (Li, M, X)
- **Basis elements** and **concentration constraints** for CE setup
- **Supercell expansion matrix** passed to `CECrystal`

---

## Script Descriptions

### `GenRandom.py` — Single-Generation Mixed Substitution

Generates structures with **simultaneous mixing on both cation (M) and anion (X) sites** at a user-specified halogen concentration `C2`. Structures are written to a fresh `.db` file (old file is deleted). Uses `generation_number=0`.

**Key parameter:** `concs=[C1, C2]` — concentration list where `C2` controls halogen substitution fraction.

---

### `GenRandomAll.py` — Multi-Generation Mixed Substitution

Functionally identical to `GenRandom.py` but designed for **iterative generation campaigns**. Accepts a `num_gen` parameter so successive generations can be appended to the same database without deletion.

**Key parameter:** `num_gen` — generation index (default `0`).

---

### `GenRandomMetal.py` — Metal-Site Substitution Only

Generates structures where **only the cation (M) site** is substituted between two metals, while the halogen site carries a fixed concentration `C1`. Uses `generation_number=0` and reinitializes the database each run.

**Key parameter:** `concs=[C1]` — single concentration controlling metal mixing.

---

### `GenRandomDefect.py` — Defect Structures

Extends the mixed substitution model to explicitly include **Li vacancies** (modeled as `'X'` placeholder on Li sites) and **anti-site defects** (M atoms on normally vacant/Li sites). The `conc` parameter (scalar) controls the defect concentration.

**Key parameter:** `conc` — scalar defect concentration (default `1`).

---

### `GenRandomDefectAll.py` — Multi-Generation Defect Structures

Multi-generation variant of `GenRandomDefect.py`. Supports the same defect model across all six templates with a `num_gen` parameter for iterative generation.

---

## Installation & Dependencies

### Prerequisites

# Python 3.8+
pip install ase
pip install clease
pip install numpy
pip install matplotlib

| Package | Purpose |
|---|---|
| `ase` | Atomic Simulation Environment — crystal structure I/O, atom manipulation |
| `clease` | Cluster Expansion for ASE — CE settings, structure generation |
| `numpy` | Numerical arrays and matrix operations |
| `matplotlib` | Plotting (imported but used for post-processing) |

### Clone the Repository

git clone https://github.com/chiku-parida/LMX_data.git
cd LMX_data

### Set Up CIF File Paths

All scripts rely on reference CIF files stored at a local path. Before running, update the `read(...)` paths in each script to point to your CIF directory:

# Example — update this line in each method:
f11 = read('/path/to/your/cif_files/LiInF.cif')

The required CIF files are:

| Template | CIF File |
|---|---|
| F11 | `LiInF.cif` |
| Cl12 | `LiScCl.cif` |
| Br12 | `LiErBr.cif` |
| F15 | `LiGaF.cif` |
| Cl62 | `LiYbCl.cif` |
| Cl164 | `LiErCl.cif` |
| F165 | `LiScF.cif` |

---

## Usage

### Basic Example — Generate Cl12 Structures (Mixed Metal + Halogen)

from GenRandom import CrysGen

gen = CrysGen(
    dbname   = "LiScCl_mixed.db",
    template = "Cl12",
    metals   = ["Sc", "Y"],       # M1, M2 — two metals to mix on cation site
    halogens = ["Cl", "Br"],      # X1, X2 — two halogens to mix on anion site
    concs    = [0.5, 0.3],        # [C1, C2] — concentration fractions
    struct_per_gen = 50,
    mixing_metals   = True,
    mixing_halogens = True,
)

gen.Cl12()


### Defect Structure Example

from GenRandomDefect import CrysGen

gen = CrysGen(
    dbname   = "LiErCl_defect.db",
    template = "Cl164",
    metals   = ["Er", "Y"],
    halogens = ["Cl", "Br"],
    conc     = 0.25,             
    num_gen  = 0,
    struct_per_gen = 20,
)

gen.Cl164()

---

## CrysGen Class 

All five scripts expose a `CrysGen` class. Below is the consolidated parameter reference.

### Constructor Parameters

| Parameter | Type | Required | Default | Description |
|---|---|---|---|---|
| `dbname` | `str` | Yes | — | Path/name of the ASE SQLite database file to write structures into |
| `template` | `str` | Yes | — | Crystal structure template key (e.g., `'Cl12'`, `'F165'`) |
| `metals` | `list[str]` | Yes | — | Two metal element symbols for the cation sublattice, e.g. `['Sc', 'Y']` |
| `halogens` | `list[str]` | Yes | — | Two halogen element symbols for the anion sublattice, e.g. `['Cl', 'F']` |
| `concs` | `list[float]` | GenRandom / GenRandomMetal | — | Concentration list `[C1, C2]`; `C2` controls halogen fraction |
| `conc` | `float` | GenRandomDefect | `1` | Scalar defect concentration for vacancy/anti-site models |
| `num_gen` | `int` | GenRandomAll / Defect variants | `0` | Generation index; allows appending to existing database |
| `struct_per_gen` | `int` | No | `100` / `5` | Number of structures to generate per call |
| `mixing_metals` | `bool` | No | `True` | Enable metal mixing on cation sites (must be `True`) |
| `mixing_halogens` | `bool` | No | `True` | Enable halogen mixing on anion sites (must be `True`) |

### Methods

Each method corresponds to a template. Call the method matching your chosen `template` string:

| Method | Template String | Space Group |
|---|---|---|
| `gen.F11()` | `'F11'` | 11 |
| `gen.Cl12()` | `'Cl12'` | 12 (Cl) |
| `gen.Br12()` | `'Br12'` | 12 (Br) |
| `gen.F15()` | `'F15'` | 15 |
| `gen.Cl62()` | `'Cl62'` | 62 |
| `gen.Cl164()` | `'Cl164'` | 164 |
| `gen.F165()` | `'F165'` | 165 |

> **Important:** The `template` string passed to the constructor must exactly match the method name. Passing a mismatched template raises a `ValueError`.

---

## Template Reference Table

| Template | SG# | Z | Li sites | M sites | X sites | Supercell |
|---|---|---|---|---|---|---|
| F11 | 11 | 2 | 4f × 4 | 2e × 2 | 2e/4f × 8 | 2×1×1 |
| Cl12 | 12 | 2 | 4g/4h/2d | 2a | 4i/8j × 2 | 2×2×2 |
| Br12 | 12 | 2 | 4g/4h | 2a | 4i/8j × 2 | 2×2×2 |
| F15 | 15 | 4 | 8f × 5 | 4a/8f | 8f × 9 | 2×1×1 |
| Cl62 | 62 | 4 | 8d × 2 | 4c | 8d/4c × 4 | 1×1×2 |
| Cl164 | 164 | 3 | 6g/6h | 1a/2d × 3 | 6i × 3 | 1×1×2 |
| F165 | 165 | 6 | 6f/12g | 2b/4d | 12g × 3 | 2×2×1 |

---

## Notes & Known Issues

- The CIF file path must be updated before use.
- **Database deletion in `GenRandom` and `GenRandomMetal`:** These scripts call `os.system(f"rm {db_name}")` at the start of each run, deleting any existing database. Use `GenRandomAll` or `GenRandomDefectAll` if you want to accumulate structures across runs.
- **Both flags must be `True`:** The `mixing_metals` and `mixing_halogens` parameters currently only support `True`. Passing `False` raises a `ValueError`.
- **CLEASE version compatibility:** The scripts use `CECrystal`, `NewStructures`, and `Concentration` from `clease`. Ensure your CLEASE version 1.0.6.

---

## Author

**CHIKU PARIDA** — [@chiku-parida](https://github.com/chiku-parida)

---

## Citation

If you use this codebase in your research, please cite the repository:

@misc{parida2025lmxdata,
  author = {Parida, Chiku},
  title  = {LMX\_data: A Database for Design and Optimization of Halide Based Solid Electrolytes for LIBs},
  year   = {2025},
  url    = {https://github.com/chiku-parida/LMX_data}
}
