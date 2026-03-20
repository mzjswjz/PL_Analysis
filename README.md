# PL Analysis

Analysis tools for photoluminescence (PL) measurements on photovoltaic and optoelectronic devices.

## Overview

Photoluminescence (PL) analysis is a non-contact optical characterization technique used to study the optical properties and quality of semiconductor materials and devices. This repository contains:

- Steady-state PL analysis
- Time-resolved PL (TRPL) analysis
- PL mapping and imaging tools
- Temperature-dependent PL analysis

## What is Photoluminescence?

When a material absorbs photons with energy greater than its bandgap, electron-hole pairs are generated. As these carriers recombine, they can emit light (photoluminescence). The PL signal provides information about:

- Bandgap energy
- Carrier lifetime
- Defect density
- Material quality
- Impurity concentrations

## Types of Analysis

### Steady-State PL
- Emission spectra at fixed excitation
- Power-dependent PL (fluence studies)
- Temperature-dependent PL

### Time-Resolved PL (TRPL)
- Carrier lifetime measurement
- Decay dynamics
- Delayed luminescence

### PL Mapping
- Spatial distribution of PL intensity
- Defect imaging
- Uniformity assessment

## Key Equations

### Radiative Recombination Rate
```
B × n × p
```

### PL Intensity
```
I_PL ∝ (E_gap) × α(hν) × φ_exc
```

### Carrier Lifetime (TRPL)
```
I(t) = I_0 × exp(-t/τ)
```

## Data Format

Standard input formats:
- CSV (wavelength, intensity)
- TXT (time, counts)
- HDF5 (multi-dimensional datasets)

## Usage

```python
python pl_analyze.py --input spectrum.csv --type steady
python trpl_analyze.py --input decay.dat --fit bi-exponential
```

## Dependencies

- NumPy
- SciPy
- Matplotlib
- pandas
- lmfit (for curve fitting)

## Author

mzjswjz
