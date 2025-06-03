
mwind in ICON: https://gitlab.dkrz.de/icon/icon-mpim/-/blob/master/src/atm_phy_aes/tmx/mo_vdf_diag_smag.f90#L162

# ICON Surface Exchange Coefficients – Winkler et al. (2025)

This repository contains the code and routines of the parameterization of surface exchange coefficients in the ICON atmospheric model. The core function implemented here is a Python version of the `sfc_exchange_coefficients` routine used in ICON.

---

## Environment Setup

We recommend using **Micromamba** for managing dependencies. To create and activate the environment:

```bash
# Create the environment from the environment.yml file
micromamba create --file environment.yml

# Activate the environment (replace 'icon_env' with your environment name if different)
micromamba activate icon_env
```

> The environment name is already specified in `environment.yml`.

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
- Diagnostic quantities like the Richardson number (RIB) and stability corrections

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

## Contact

For questions or feedback, please open an issue [here in this repository](https://github.com/mariuswinkler/Winkler_et_al_ICON_Surface_Exchange_paper_2025/issues).
