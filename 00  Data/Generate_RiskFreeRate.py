import pandas as pd
import numpy  as np
from pathlib import Path


# ============================================================
# 0. 路徑設定
# ============================================================

data_dir = Path(
    r"D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel\Data"
)

output_dir = Path(
    r"D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel\Code\00  Output"
)

input_file = data_dir / "raw_ZeroCouponYield_19962025.csv"
output_file = output_dir / "RiskFreeRate19962025.txt"


# ============================================================
# 1. 讀取 csv
# ============================================================

df = pd.read_csv(input_file)

# 欄位名稱保險處理：去掉前後空白、轉小寫
df.columns = df.columns.str.strip().str.lower()

print("CSV columns:", df.columns.tolist())


# ============================================================
# 2. 自動判斷欄位名稱
# ============================================================
# 你的資料通常會是：
# date / days / rate
# 或 caldt / days / rate
# 若欄位名稱不同，可以直接在這裡手動指定。

date_candidates = ["date", "caldt", "datadate"]
days_candidates = ["days", "ttm", "maturity", "term"]
rate_candidates = ["rate", "yield", "zero_rate", "zerorate", "yld"]

date_col = next((c for c in date_candidates if c in df.columns), None)
days_col = next((c for c in days_candidates if c in df.columns), None)
rate_col = next((c for c in rate_candidates if c in df.columns), None)

if date_col is None:
    raise ValueError(f"找不到日期欄位，請檢查欄位名稱：{df.columns.tolist()}")

if days_col is None:
    raise ValueError(f"找不到天數欄位，請檢查欄位名稱：{df.columns.tolist()}")

if rate_col is None:
    raise ValueError(f"找不到利率欄位，請檢查欄位名稱：{df.columns.tolist()}")

print(f"date column = {date_col}")
print(f"days column = {days_col}")
print(f"rate column = {rate_col}")


# ============================================================
# 3. 整理日期、天數、利率
# ============================================================

out = df[[date_col, days_col, rate_col]].copy()
out.columns = ["date", "days", "rate"]

# 日期轉成 YYYYMMDD
# 可處理 19960104、1996-01-04、datetime 等格式
out["date"] = pd.to_datetime(out["date"].astype(str), errors="coerce")
out["date"] = out["date"].dt.strftime("%Y%m%d")

# days / rate 轉數值
out["days"] = pd.to_numeric(out["days"], errors="coerce")
out["rate"] = pd.to_numeric(out["rate"], errors="coerce")

# 移除缺值
out = out.dropna(subset=["date", "days", "rate"])

# days 轉整數
out["days"] = out["days"].astype(int)

# rate 除以 100
out["rate"] = out["rate"] / 100

# 日期小到大；同一天內 days 小到大
out = out.sort_values(["date", "days"]).reset_index(drop=True)


# ============================================================
# 4. 輸出成 txt：無標題列、無 index、空白分隔
# ============================================================

output_dir.mkdir(parents=True, exist_ok=True)

with open(output_file, "w", encoding="utf-8", newline="\n") as f:
    for _, row in out.iterrows():
        date = row["date"]
        days = int(row["days"])
        rate = format(row["rate"], ".10g")

        f.write(f"{date} {days} {rate}\n")


# ============================================================
# 5. 檢查輸出
# ============================================================

print(f"Saved: {output_file}")
print(f"Number of rows: {len(out):,}")

print("\nHead:")
print(out.head(10))

print("\nTail:")
print(out.tail(10))


# ============================================================
# 6. 生成各 TTM 對應的 risk-free gross factor table
# ============================================================
#
# 前面已經有：
# out = DataFrame with columns ["date", "days", "rate"]
#
# date: YYYYMMDD
# days: zero-coupon maturity days
# rate: annualized decimal rate, e.g. 0.052
#
# 這裡會輸出：
# 1. Risk_Free_GrossFactor_ByTargetTTM.csv
# 2. Risk_Free_GrossFactor_ByTargetTTM_Long.csv
# 3. Risk_Free_GrossFactor_ByTargetTTM_Audit.csv
# 4. Risk_Free_GrossFactor_ByTargetTTM_Summary.csv


# ============================================================
# 6.1 設定
# ============================================================

DAY_COUNT = 365

TTM_LIST = [30, 60, 90, 180]

# 跟新版 MATLAB 主程式一致：
# rf_days = exdate - date - 1
RF_DAYS_MODE = "matlab_actual"

# 如果之後想完全仿照舊版固定天數，可以改成：
# RF_DAYS_MODE = "legacy_fixed"
LEGACY_FIXED_RF_DAYS = {
    30: 29,
    60: 59,
    90: 89,
    180: 179,
}

wide_output_file    = output_dir / "Risk_Free_GrossFactor_ByTargetTTM.csv"
long_output_file    = output_dir / "Risk_Free_GrossFactor_ByTargetTTM_Long.csv"
audit_output_file   = output_dir / "Risk_Free_GrossFactor_ByTargetTTM_Audit.csv"
summary_output_file = output_dir / "Risk_Free_GrossFactor_ByTargetTTM_Summary.csv"


# ============================================================
# 6.2 整理 RF curve table
# ============================================================

rf_curve = out.copy()

rf_curve["date"] = (
    rf_curve["date"]
    .astype(str)
    .str.replace(r"\.0$", "", regex=True)
    .str.zfill(8)
)

rf_curve["days"] = pd.to_numeric(rf_curve["days"], errors="coerce")
rf_curve["rate"] = pd.to_numeric(rf_curve["rate"], errors="coerce")

rf_curve = rf_curve.dropna(subset=["date", "days", "rate"])
rf_curve["days"] = rf_curve["days"].astype(int)

# 同一天同一個 days 若有重複，取平均
# 對齊 MATLAB:
# [Data_TTM_U, ~, idx_group] = unique(Data_TTM);
# Data_RF_U = accumarray(idx_group, Data_RF, [], @mean);
rf_curve = (
    rf_curve
    .groupby(["date", "days"], as_index=False)["rate"]
    .mean()
    .sort_values(["date", "days"])
    .reset_index(drop=True)
)

rf_curve_by_date = {}

for date, sub in rf_curve.groupby("date"):

    sub = sub.sort_values("days")

    rf_curve_by_date[date] = {
        "days": sub["days"].to_numpy(dtype=float),
        "rate": sub["rate"].to_numpy(dtype=float),
    }
    
    
# ============================================================
# 6.3 Interpolation / extrapolation functions
# ============================================================

def linear_interp_extrap(x, y, x0):
    """
    Replicate MATLAB interp1(x, y, x0, 'linear', 'extrap').

    x: available zero-coupon maturity days
    y: annualized risk-free rates
    x0: target maturity days
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x0 = float(x0)

    valid = np.isfinite(x) & np.isfinite(y)

    x = x[valid]
    y = y[valid]

    if len(x) == 0:
        return np.nan, "missing_curve"

    if len(x) == 1:
        return float(y[0]), "single_point"

    idx_sort = np.argsort(x)

    x = x[idx_sort]
    y = y[idx_sort]

    if x0 < x[0]:

        slope = (y[1] - y[0]) / (x[1] - x[0])
        value = y[0] + slope * (x0 - x[0])

        return float(value), "extrap_below"

    if x0 > x[-1]:

        slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
        value = y[-1] + slope * (x0 - x[-1])

        return float(value), "extrap_above"

    value = np.interp(x0, x, y)

    return float(value), "interpolate"


def get_rf_for_target_date(target_date, rf_days, max_lag_days=10):
    """
    Replicate Function__RF_TTM logic:

    1. Search target date first.
    2. If not found, search up to 10 calendar days before.
    3. Interpolate annualized RF by target rf_days.
    4. If rf_days is outside available maturity range, linearly extrapolate.
    """

    target_date_str = str(int(target_date)).zfill(8)
    target_dt = pd.to_datetime(target_date_str, format="%Y%m%d")

    for lag in range(max_lag_days + 1):

        search_dt = target_dt - pd.Timedelta(days=lag)
        search_date_str = search_dt.strftime("%Y%m%d")

        if search_date_str in rf_curve_by_date:

            curve = rf_curve_by_date[search_date_str]

            rf_annualized, interp_method = linear_interp_extrap(
                curve["days"],
                curve["rate"],
                rf_days,
            )

            return {
                "rf_annualized": rf_annualized,
                "rf_source_date": int(search_date_str),
                "rf_lag_days": lag,
                "interp_method": interp_method,
                "status": "Success" if np.isfinite(rf_annualized) else "Missing RF",
            }

    return {
        "rf_annualized": np.nan,
        "rf_source_date": np.nan,
        "rf_lag_days": np.nan,
        "interp_method": "missing_rf_curve",
        "status": "Missing RF curve within 10 days",
    }


# ============================================================
# 6.4 For each TTM, compute RF gross factor
# ============================================================

wide_tables = []
long_tables = []
audit_tables = []
summary_rows = []

for target_ttm in TTM_LIST:

    ttm_file = output_dir / f"TTM_{target_ttm}.csv"

    if not ttm_file.exists():

        print(f"[SKIP] Missing TTM file: {ttm_file}")
        continue

    target = pd.read_csv(ttm_file)
    target.columns = target.columns.str.strip().str.lower()

    required_cols = {"date", "exdate"}
    missing_cols = required_cols - set(target.columns)

    if missing_cols:
        raise ValueError(f"{ttm_file.name} is missing required columns: {missing_cols}")

    target["date"] = pd.to_numeric(target["date"], errors="coerce")
    target["exdate"] = pd.to_numeric(target["exdate"], errors="coerce")

    target = target.dropna(subset=["date", "exdate"])

    target["date"] = target["date"].astype(int)
    target["exdate"] = target["exdate"].astype(int)

    target = target.sort_values("date").reset_index(drop=True)

    date_dt = pd.to_datetime(target["date"].astype(str), format="%Y%m%d")
    exdate_dt = pd.to_datetime(target["exdate"].astype(str), format="%Y%m%d")

    target["calendar_ttm"] = (exdate_dt - date_dt).dt.days

    if RF_DAYS_MODE == "matlab_actual":

        # 跟新版 MATLAB RND 主程式一致：
        # Data(:, Index_TTM) = Exp_Corrected_Num_Selected - Date_Num_Selected;
        # Data(:, Index_TTM) = Data(:, Index_TTM) - 1;
        target["rf_days"] = target["calendar_ttm"] - 1

        # 欄位名稱不要寫死 29d，因為 Juneteenth 等特殊情況可能是 30d
        rf_col = f"rf_gross_TTM{target_ttm}"

    elif RF_DAYS_MODE == "legacy_fixed":

        target["rf_days"] = LEGACY_FIXED_RF_DAYS[target_ttm]
        rf_col = f"rf_gross_{LEGACY_FIXED_RF_DAYS[target_ttm]}d"

    else:

        raise ValueError("RF_DAYS_MODE must be either 'matlab_actual' or 'legacy_fixed'.")

    if (target["rf_days"] <= 0).any():

        bad = target[target["rf_days"] <= 0]

        raise ValueError(
            f"Non-positive rf_days found in TTM_{target_ttm}.csv:\n{bad}"
        )

    records = []

    for _, row in target.iterrows():

        rf_info = get_rf_for_target_date(
            target_date=row["date"],
            rf_days=row["rf_days"],
            max_lag_days=10,
        )

        rf_annualized = rf_info["rf_annualized"]

        if np.isfinite(rf_annualized):
            rf_gross = np.exp(rf_annualized * row["rf_days"] / DAY_COUNT)
        else:
            rf_gross = np.nan

        records.append({
            "target_ttm": target_ttm,
            "date": int(row["date"]),
            "exdate": int(row["exdate"]),
            "calendar_ttm": int(row["calendar_ttm"]),
            "rf_days": int(row["rf_days"]),
            "rf_annualized": rf_annualized,
            "rf_gross": rf_gross,
            "rf_source_date": rf_info["rf_source_date"],
            "rf_lag_days": rf_info["rf_lag_days"],
            "interp_method": rf_info["interp_method"],
            "status": rf_info["status"],
        })

    result = pd.DataFrame(records)

    valid = result[result["status"] == "Success"].copy()

    # Wide format：仿照舊 Risk_Free_Rate.csv
    wide_piece = pd.DataFrame({
        f"date_{target_ttm}": valid["date"].astype("Int64").reset_index(drop=True),
        rf_col: valid["rf_gross"].reset_index(drop=True),
    })

    wide_tables.append(wide_piece)

    # Long format：比較適合 merge / debug
    long_piece = valid[[
        "target_ttm",
        "date",
        "exdate",
        "calendar_ttm",
        "rf_days",
        "rf_annualized",
        "rf_gross",
    ]].copy()

    long_tables.append(long_piece)

    audit_tables.append(result)

    summary_rows.append({
        "target_ttm": target_ttm,
        "n_target_rows": len(result),
        "n_success": int((result["status"] == "Success").sum()),
        "n_missing": int((result["status"] != "Success").sum()),
        "min_date": int(result["date"].min()) if len(result) > 0 else np.nan,
        "max_date": int(result["date"].max()) if len(result) > 0 else np.nan,
        "min_rf_days": int(result["rf_days"].min()) if len(result) > 0 else np.nan,
        "max_rf_days": int(result["rf_days"].max()) if len(result) > 0 else np.nan,
        "n_interpolate": int((result["interp_method"] == "interpolate").sum()),
        "n_extrap_below": int((result["interp_method"] == "extrap_below").sum()),
        "n_extrap_above": int((result["interp_method"] == "extrap_above").sum()),
        "n_lagged_rf_date": int((result["rf_lag_days"] > 0).sum()),
    })

    print(f"TTM = {target_ttm}")
    print(f"  target rows        = {len(result):,}")
    print(f"  success rows       = {len(valid):,}")
    print(f"  missing rows       = {(result['status'] != 'Success').sum():,}")
    print(f"  rf_days range      = {result['rf_days'].min()} to {result['rf_days'].max()}")
    print(f"  lagged RF rows     = {(result['rf_lag_days'] > 0).sum():,}")
    print(f"  output column name = {rf_col}")
    print()


# ============================================================
# 6.5 Save outputs
# ============================================================

if wide_tables:

    df_rf_wide = pd.concat(wide_tables, axis=1)
    df_rf_wide.to_csv(wide_output_file, index=False, float_format="%.15g")

else:

    df_rf_wide = pd.DataFrame()


if long_tables:

    df_rf_long = pd.concat(long_tables, axis=0, ignore_index=True)
    df_rf_long.to_csv(long_output_file, index=False, float_format="%.15g")

else:

    df_rf_long = pd.DataFrame()


if audit_tables:

    df_rf_audit = pd.concat(audit_tables, axis=0, ignore_index=True)
    df_rf_audit.to_csv(audit_output_file, index=False, float_format="%.15g")

else:

    df_rf_audit = pd.DataFrame()


df_rf_summary = pd.DataFrame(summary_rows)
df_rf_summary.to_csv(summary_output_file, index=False)


# ============================================================
# 6.6 Check outputs
# ============================================================

print("Saved files:")
print(f"  Wide format : {wide_output_file}")
print(f"  Long format : {long_output_file}")
print(f"  Audit file  : {audit_output_file}")
print(f"  Summary     : {summary_output_file}")

print("\nSummary:")
print(df_rf_summary)

print("\nWide format head:")
print(df_rf_wide.head())

print("\nWide format tail:")
print(df_rf_wide.tail())