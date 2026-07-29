"""Assemble g(R), g'(R), and g''(R) from the compact CSV files.

Required packages: numpy and pandas.

The derivatives are taken with respect to the gross return R. For a selected
TTM, rows of the resulting matrices correspond to quote dates and columns
correspond to the common 30,000-point gross-return grid.
"""

from pathlib import Path

import numpy as np
import pandas as pd


# Replace this value with the downloaded "19  Output - g" folder.
DATA_DIR = Path(r"REPLACE_WITH_DOWNLOADED_19_OUTPUT_G_FOLDER")

# Change this value to 60, 90, or 180 when needed.
TTM = 30


def load_compact_components(data_dir, ttm):
    """Load and validate the compact components for one TTM."""
    r_grid = pd.read_csv(data_dir / "CubicBSpline_R_Grid.csv")
    shape_all = pd.read_csv(
        data_dir / "CubicBSpline_GShape_By_TTM.csv"
    )
    time_scale_all = pd.read_csv(
        data_dir / "CubicBSpline_GTimeScale_By_TTM.csv"
    )

    shape = (
        shape_all.loc[shape_all["TTM"] == ttm]
        .sort_values("grid_index")
        .reset_index(drop=True)
    )
    time_scale = (
        time_scale_all.loc[time_scale_all["TTM"] == ttm]
        .sort_values("quote_date")
        .reset_index(drop=True)
    )

    if shape.empty:
        raise ValueError(f"No g-shape data were found for TTM {ttm}.")
    if time_scale.empty:
        raise ValueError(f"No time-scale data were found for TTM {ttm}.")
    if len(r_grid) != len(shape):
        raise ValueError("The R grid and g-shape have different row counts.")
    if not np.array_equal(r_grid["grid_index"], shape["grid_index"]):
        raise ValueError("The R grid and g-shape grid indexes do not match.")
    if time_scale["quote_date"].duplicated().any():
        raise ValueError("Each quote date must have exactly one time scale.")

    return r_grid, shape, time_scale


def assemble_all_months(r_grid, shape, time_scale):
    """Assemble all quote dates as three date-by-return matrices."""
    scale = time_scale["g_scale"].to_numpy(dtype=float)[:, None]
    g = scale * shape["g_shape"].to_numpy(dtype=float)[None, :]
    g_d1 = scale * shape["g_shape_d1"].to_numpy(dtype=float)[None, :]
    g_d2 = scale * shape["g_shape_d2"].to_numpy(dtype=float)[None, :]

    quote_dates = time_scale["quote_date"].to_numpy(dtype=np.int64)
    expiration_dates = time_scale["expiration_date"].to_numpy(dtype=np.int64)
    gross_return = r_grid["gross_return_R"].to_numpy(dtype=float)

    return {
        "quote_dates": quote_dates,
        "expiration_dates": expiration_dates,
        "gross_return_R": gross_return,
        "g": g,
        "g_d1": g_d1,
        "g_d2": g_d2,
    }


def assemble_one_month(r_grid, shape, time_scale, quote_date):
    """Assemble g and its derivatives for one quote date."""
    selected = time_scale.loc[time_scale["quote_date"] == quote_date]
    if len(selected) != 1:
        raise ValueError(f"Expected one time-scale row for {quote_date}.")

    scale = float(selected["g_scale"].iloc[0])
    return pd.DataFrame(
        {
            "grid_index": r_grid["grid_index"].to_numpy(dtype=np.int32),
            "gross_return_R": r_grid["gross_return_R"].to_numpy(dtype=float),
            "g": scale * shape["g_shape"].to_numpy(dtype=float),
            "g_d1": scale * shape["g_shape_d1"].to_numpy(dtype=float),
            "g_d2": scale * shape["g_shape_d2"].to_numpy(dtype=float),
        }
    )


r_grid, shape, time_scale = load_compact_components(DATA_DIR, TTM)
all_months = assemble_all_months(r_grid, shape, time_scale)

print(f"TTM: {TTM}")
print(f"Number of quote dates: {len(all_months['quote_dates'])}")
print(f"Number of R-grid points: {len(all_months['gross_return_R'])}")
print(f"g matrix shape: {all_months['g'].shape}")
print(f"g_d1 matrix shape: {all_months['g_d1'].shape}")
print(f"g_d2 matrix shape: {all_months['g_d2'].shape}")

# Example: extract one quote date as a long-form DataFrame.
example_quote_date = int(all_months["quote_dates"][0])
one_month = assemble_one_month(
    r_grid,
    shape,
    time_scale,
    example_quote_date,
)

print(f"\nExample quote date: {example_quote_date}")
print(one_month.head())
