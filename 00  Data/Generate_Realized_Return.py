import pandas as pd
from pathlib import Path


# =============================================================================
# Paths
# =============================================================================

data_dir = Path(
    r"D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel\Data"
)

output_dir = Path(
    r"D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel\Code\00  Output"
)

spx_file = data_dir / "raw_SPX_19962025.csv"

ttm_list = [30, 60, 90, 180]


# =============================================================================
# Load SPX closing price
# =============================================================================
# Expected raw SPX format:
# secid,date,close
# 108105,19960102,620.73
# ...

spx = pd.read_csv(spx_file)
spx.columns = spx.columns.str.strip().str.lower()

required_cols = {"secid", "date", "close"}
missing_cols = required_cols - set(spx.columns)

if missing_cols:
    raise ValueError(f"SPX file is missing required columns: {missing_cols}")

# Keep SPX index only
spx = spx[spx["secid"] == 108105].copy()

spx["date"] = pd.to_numeric(spx["date"], errors="coerce")
spx["close"] = pd.to_numeric(spx["close"], errors="coerce")

spx = spx.dropna(subset=["date", "close"])
spx["date"] = spx["date"].astype(int)

spx = spx.sort_values("date").drop_duplicates(subset=["date"], keep="last")

price_map = spx.set_index("date")["close"]

min_price_date = int(spx["date"].min())
max_price_date = int(spx["date"].max())

print("SPX price range:")
print(f"  min date = {min_price_date}")
print(f"  max date = {max_price_date}")
print(f"  number of trading days = {len(spx):,}")
print()


# =============================================================================
# Generate realized returns for each TTM file
# =============================================================================

summary_rows = []

for ttm in ttm_list:

    ttm_file = output_dir / f"TTM_{ttm}.csv"
    realized_file = output_dir / f"Realized_Return_TTM_{ttm}.csv"
    audit_file = output_dir / f"Realized_Return_TTM_{ttm}_audit.csv"

    if not ttm_file.exists():
        print(f"[SKIP] TTM file not found: {ttm_file}")
        continue

    target = pd.read_csv(ttm_file)
    target.columns = target.columns.str.strip().str.lower()

    required_ttm_cols = {"date", "exdate"}
    missing_ttm_cols = required_ttm_cols - set(target.columns)

    if missing_ttm_cols:
        raise ValueError(f"{ttm_file.name} is missing required columns: {missing_ttm_cols}")

    target["date"] = pd.to_numeric(target["date"], errors="coerce")
    target["exdate"] = pd.to_numeric(target["exdate"], errors="coerce")

    target = target.dropna(subset=["date", "exdate"])
    target["date"] = target["date"].astype(int)
    target["exdate"] = target["exdate"].astype(int)

    target = target.sort_values("date").reset_index(drop=True)

    # Attach quote-date and exdate closing prices
    target["price_date"] = target["date"].map(price_map)
    target["price_exdate"] = target["exdate"].map(price_map)

    # Diagnostics flags
    target["has_price_date"] = target["price_date"].notna()
    target["has_price_exdate"] = target["price_exdate"].notna()

    target["drop_reason"] = ""

    target.loc[~target["has_price_date"], "drop_reason"] = "missing quote-date SPX close"
    target.loc[~target["has_price_exdate"], "drop_reason"] = "missing exdate SPX close"

    target.loc[
        (~target["has_price_date"]) & (~target["has_price_exdate"]),
        "drop_reason"
    ] = "missing both quote-date and exdate SPX close"

    # Main valid sample
    valid = target[target["has_price_date"] & target["has_price_exdate"]].copy()

    valid["realized_ret"] = valid["price_exdate"] / valid["price_date"]

    realized = valid[["date", "realized_ret"]].copy()

    # Save realized return
    realized.to_csv(realized_file, index=False, float_format="%.15g")

    # Save audit file
    audit_cols = [
        "date",
        "exdate",
        "price_date",
        "price_exdate",
        "has_price_date",
        "has_price_exdate",
        "drop_reason",
    ]

    target[audit_cols].to_csv(audit_file, index=False, float_format="%.15g")

    n_total = len(target)
    n_valid = len(realized)
    n_drop = n_total - n_valid

    first_valid_date = int(realized["date"].min()) if n_valid > 0 else None
    last_valid_date = int(realized["date"].max()) if n_valid > 0 else None

    first_dropped_date = None
    if n_drop > 0:
        dropped = target[target["drop_reason"] != ""].copy()
        first_dropped_date = int(dropped["date"].min())

    summary_rows.append({
        "TTM": ttm,
        "n_target_rows": n_total,
        "n_realized_rows": n_valid,
        "n_dropped_rows": n_drop,
        "first_valid_date": first_valid_date,
        "last_valid_date": last_valid_date,
        "first_dropped_date": first_dropped_date,
        "output_file": realized_file.name,
        "audit_file": audit_file.name,
    })

    print(f"TTM = {ttm}")
    print(f"  target rows   = {n_total:,}")
    print(f"  realized rows = {n_valid:,}")
    print(f"  dropped rows  = {n_drop:,}")
    print(f"  saved         = {realized_file}")
    print(f"  audit saved   = {audit_file}")

    if n_drop > 0:
        print("  dropped observations:")
        print(target.loc[target["drop_reason"] != "", ["date", "exdate", "drop_reason"]].tail(10))

    print()


# =============================================================================
# Save summary
# =============================================================================

summary = pd.DataFrame(summary_rows)

summary_file = output_dir / "Realized_Return_Summary.csv"
summary.to_csv(summary_file, index=False)

print("Summary:")
print(summary)
print()
print(f"Saved summary: {summary_file}")