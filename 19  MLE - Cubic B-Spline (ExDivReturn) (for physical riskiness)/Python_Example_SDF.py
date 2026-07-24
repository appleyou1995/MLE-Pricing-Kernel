"""
Load or reconstruct the estimated cubic B-spline SDF in Python.
Required packages: numpy, pandas, and scipy.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import BSpline


# User setting: path to the downloaded output folder.
DATA_DIR = Path(r"REPLACE_WITH_DOWNLOADED_FOLDER")

# Example setting. Change TTM to 60, 90, or 180 when needed.
TTM = 30
METHOD = "reconstruct"  # Use "reconstruct" or "read".

PREFIX = f"CubicBSpline_TTM{TTM:03d}"
DATE_COLUMN = f"date_TTM_{TTM}"
EXDATE_COLUMN = f"exdate_TTM_{TTM}"
KAPPA_COLUMN = f"kappa_TTM_{TTM}"


# Method 1: Reconstruct the SDF from theta, knots, and monthly kappa.
theta = pd.read_csv(DATA_DIR / f"{PREFIX}_Theta.csv")["theta"].to_numpy(float)

knots = (
    pd.read_csv(DATA_DIR / f"{PREFIX}_Knots.csv")
    .sort_values("knot_index")["knot_value"]
    .to_numpy(float)
)

kappa_table = (
    pd.read_csv(DATA_DIR / f"{PREFIX}_Kappa.csv")
    .dropna(subset=[DATE_COLUMN, EXDATE_COLUMN, KAPPA_COLUMN])
    .sort_values(DATE_COLUMN)
    .reset_index(drop=True)
)

if kappa_table[DATE_COLUMN].duplicated().any():
    raise ValueError("Each quote date must have exactly one kappa value.")

degree = 3
number_of_basis_functions = len(knots) - degree - 1
if len(theta) != number_of_basis_functions:
    raise ValueError("The numbers of theta coefficients and basis functions differ.")

# This is the same 30,000-point gross-return grid used in the MATLAB output.
gross_return = np.linspace(0.003, 3.0, 30_000)
basis_matrix = BSpline.design_matrix(
    gross_return,
    knots,
    k=degree,
    extrapolate=False,
).toarray()

spline_component = basis_matrix @ theta
grid_index = np.arange(1, len(gross_return) + 1, dtype=np.int32)


def reconstruct_monthly_sdf(month_row):
    """Return the full 30,000-point SDF for one row of the kappa table."""
    quote_date = int(month_row[DATE_COLUMN])
    kappa_t = float(month_row[KAPPA_COLUMN])
    sdf_value = np.exp(kappa_t + spline_component)

    return pd.DataFrame(
        {
            "target_ttm_days": TTM,
            "quote_date": quote_date,
            "grid_index": grid_index,
            "gross_return_R": gross_return,
            "sdf_M": sdf_value,
        }
    )


expected_rows = len(kappa_table) * len(gross_return)


def reconstruct_all_sdf():
    """Reconstruct every quote date listed in the kappa file."""
    result = pd.concat(
        (reconstruct_monthly_sdf(row) for _, row in kappa_table.iterrows()),
        ignore_index=True,
    )
    if len(result) != expected_rows:
        raise ValueError(f"Expected {expected_rows:,} reconstructed SDF rows.")
    return result


# Method 2: Read all ready-made annual long-form R-SDF files directly.
def read_all_ready_made_sdf():
    """Read every annual long-form SDF file for the selected TTM."""
    annual_grid_files = sorted(
        (DATA_DIR / "SDF_Grid_CSV").glob(f"{PREFIX}_SDFGrid_*.csv")
    )
    if not annual_grid_files:
        raise FileNotFoundError(
            f"No annual SDF grid files were found for TTM {TTM}."
        )

    result = pd.concat(
        (pd.read_csv(file) for file in annual_grid_files),
        ignore_index=True,
    )
    if len(result) != expected_rows:
        raise ValueError(f"Expected {expected_rows:,} ready-made SDF rows.")
    return result


# Both methods return all quote dates in one long-form DataFrame named all_sdf.
# Running only one method at a time avoids holding two large copies in memory.
if METHOD == "reconstruct":
    all_sdf = reconstruct_all_sdf()
elif METHOD == "read":
    all_sdf = read_all_ready_made_sdf()
else:
    raise ValueError('METHOD must be either "reconstruct" or "read".')

print(f"Method: {METHOD}")
print(f"Number of quote dates: {all_sdf['quote_date'].nunique()}")
print(f"Number of rows: {len(all_sdf):,}")
print(all_sdf.head())
