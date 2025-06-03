## ICON Surface Exchange Coefficients (Winkler et al., 2025 [in prep.])

This repository contains the code and routines of the parameterization of surface exchange coefficients in the ICON atmospheric model. The core function implemented here is a Python version of the `sfc_exchange_coefficients` routine used in ICON.

---

## Environment Setup

We recommend using **Micromamba** for managing dependencies. To create and activate the environment:

```bash
# Create the environment from the environment.yml file
micromamba create --file environment.yml

# Activate the environment
micromamba activate ICON_sfc_exchange_coefficients
```

---

## Core Function: `sfc_exchange_coefficients`

Located in [`src/functions.py`](src/functions.py), this function computes surface exchange coefficients for **momentum** and **heat** based on Monin–Obukhov similarity theory.

### Inputs:
- Atmospheric state at the lowest model level (wind, potential temperature, humidity)
- Surface values (surface potential temperature, saturation humidity)
- Surface roughness and reference height

### Outputs:
- `cD`, `cH`: Exchange coefficients for momentum and heat
- `cD_neutral`, `cH_neutral`: Corresponding coefficients under neutral conditions
- Diagnostic quantities like the Richardson number (RIB) and stability functions

The function mirrors the iterative procedure found in the ICON routine `sfc_exchange_coefficients` (see `icon-mpim/src/atm_phy_aes/tmx/mo_vdf_diag_smag.f90`) to resolve surface-layer fluxes under varying stability.

---

## Repository Structure

```
├── environment.yml          # Micromamba environment file
├── README.md                # This file
├── src/
│   ├── functions.py         # Core implementation of sfc_exchange_coefficients
│   └── constants.py         # Physical and model constants
├── scripts/
│   └── compute_sfc_exchange_coefficient.py  # Main script for computation and plotting
├── figs/                    # Generated figures
└── .gitignore
```

---

## Usage

You can run the surface exchange coefficient computation via:

```bash
python scripts/compute_sfc_exchange_coefficient.py
```

This script loops over a range of wind speeds, calculates exchange coefficients, and generates corresponding plots for different wind thresholds and atmospheric conditions.

---

## Relevant Literature

- **Kitamura, Y., & Ito, J. (2016)**  
  [*Revisiting the Bulk Relation for Heat Flux in the Free Convection Limit*](https://doi.org/10.1007/s10546-015-0075-z), *Boundary-Layer Meteorology*, 158, 93–103.

- **Zeng, X., Zhang, Q., Johnson, D., & Tao, W. (2002)**  
  [*Parameterization of Wind Gustiness for the Computation of Ocean Surface Fluxes at Different Spatial Scales*](https://doi.org/10.1175/1520-0493(2002)130<2125:POWGFT>2.0.CO;2), *Monthly Weather Review*, 130, 2125–2133.

---

## Contact

For questions or feedback, please open an issue [here in this repository](https://github.com/mariuswinkler/Winkler_et_al_ICON_Surface_Exchange_paper_2025/issues).
