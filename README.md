# PL Analysis

Python tools for analyzing photoluminescence (PL) spectroscopy data.

## Files

| File | Description |
|------|-------------|
| `main.py` | Example usage script showing how to use the PL analysis class |
| `PL_functions/PL_analysis_tools.py` | `Photoluminescence` class — core analysis functionality |
| `PL_functions/__init__.py` | Module initializer |

## PL_analysis_tools.py Functions

The `Photoluminescence` class provides methods for loading and visualizing PL data:

| Method | Description |
|--------|-------------|
| `__init__(bin_file_id)` | Load PL data from a CSV index file mapping molecule names to data files |
| `plot_unnorm_PL(...)` | Plot un-normalized PL spectra |
| `plot_unnorm_smooth_PL(smooth_level, ...)` | Plot smoothed un-normalized PL spectra |
| `plot_norm_PL(...)` | Plot normalized PL spectra (max = 1) |
| `plot_norm_smooth_PL(smooth_level, ...)` | Plot smoothed normalized PL spectra |
| `plot_smooth_norm_PL(smooth_level, ...)` | Alternative smoothing + normalization plot |
| `convert_to_eV(wavelength)` | Convert wavelength (nm) to energy (eV) |
| `plot_PL(...)` | Base plotting method with peak finding capability |

## Features

- Load multiple PL spectra from a CSV index file
- Plot in wavelength (nm) or energy (eV) units
- Optional Savitzky-Golay smoothing
- Automatic peak detection with prominence filtering
- Normalized or absolute intensity plots
- Save plots with timestamps

## Input Format

The `bin_file_id` CSV file should contain:
```csv
Molecule_name,File_ID
Sample1,/path/to/sample1_data.txt
Sample2,/path/to/sample2_data.txt
```

Each data file should contain wavelength and counts columns.

## Usage

```python
from PL_functions.PL_analysis_tools import Photoluminescence

# Load data
pl = Photoluminescence('path/to/index.csv')

# Plot normalized spectra with peak detection
pl.plot_norm_PL(plot_energy=True, plot_peaks=True, savefig=True)
```

Or run the example:
```bash
python main.py
```

---
*Author: mzjswjz*
