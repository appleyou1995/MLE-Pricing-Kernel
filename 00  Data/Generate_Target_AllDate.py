import pandas as pd
import numpy as np
from pathlib import Path
from calendar import monthcalendar, FRIDAY
import os
import re
import exchange_calendars as xcals


# ============================================================
# 0. Parameters
# ============================================================

START_QUOTE_MONTH = "1996-01"
END_QUOTE_MONTH   = "2025-08"

END_BY = "quote_month"
END_EXP_MONTH = "2025-08"   # 這行在 END_BY = "quote_month" 時不會用到

OUTPUT_DIR = Path(
    r"D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel\Code\00  Output"
)

SAVE_AUDIT = True

PROJECT_ROOT = OUTPUT_DIR.parent.parent
OUTPUT_DIR = Path(os.environ.get("MLE_TARGET_OUTPUT_DIR", OUTPUT_DIR))
RAW_OPTION_DIR = (
    PROJECT_ROOT
    / "Data"
    / "IndexOptions1996202508_SP500"
    / "IV-Based"
)

# TTM label 對應到幾個月後的 traditional monthly SPX expiration
# 注意：這裡的 60/90/180 是沿用舊檔命名，不是嚴格 calendar days。
TTM_MONTH_MAP = {
    30: 1,
    60: 2,
    90: 3,
    180: 6,
}

# Keep the original calendar construction for 30/60 days.  For 90/180 days,
# choose an actually traded standard-monthly chain close to the original quote
# date.  The maturity bounds prevent a missing 180-day contract from silently
# turning into an approximately 90-day observation.
DATA_DRIVEN_TTM_LABELS = {90, 180}
MAX_QUOTE_DATE_SHIFT_DAYS = 14
QUOTE_DATE_SHIFT_PENALTY = 2
MIN_PRELIMINARY_OBS = 14
MIN_PRELIMINARY_UNIQUE_K = 10
TTM_BOUNDS = {
    90: (75, 105),
    180: (140, 220),
}


# ============================================================
# 1. NYSE trading calendar
# ============================================================

cal = xcals.get_calendar(
    "XNYS",
    start="1990-01-01",
    end="2030-12-31"
)


def to_date(date):
    return pd.Timestamp(date).normalize()


def is_trading_day(date):
    date = to_date(date)
    return cal.is_session(date)


def previous_trading_day(date):
    date = to_date(date)

    while not is_trading_day(date):
        date -= pd.Timedelta(days=1)

    return date


# ============================================================
# 2. Calendar utilities
# ============================================================

def third_friday(year, month):
    cal_month = monthcalendar(year, month)
    fridays = [week[FRIDAY] for week in cal_month if week[FRIDAY] != 0]

    return pd.Timestamp(year=year, month=month, day=fridays[2])


def monthly_spx_settlement_date(period_month):
    """
    Traditional monthly SPX settlement date:
    usually third Friday, but shifted backward if third Friday is not a trading day.
    """
    third_fri = third_friday(period_month.year, period_month.month)
    exdate = previous_trading_day(third_fri)

    return exdate, third_fri


def target_quote_date_from_quote_month(quote_month):
    """
    沿用舊版 Target_AllDate 的邏輯：

    quote month 的 date
    = quote month 下一個 monthly settlement date - 30 calendar days
    = 若不是交易日，往前移到前一個 NYSE trading day
    """
    base_exp_month = quote_month + 1

    base_exdate, base_third_friday = monthly_spx_settlement_date(base_exp_month)

    raw_quote_date = base_exdate - pd.Timedelta(days=30)
    quote_date = previous_trading_day(raw_quote_date)

    return {
        "quote_date": quote_date,
        "raw_quote_date": raw_quote_date,
        "base_exp_month": base_exp_month,
        "base_exdate": base_exdate,
        "base_third_friday": base_third_friday,
        "base_exdate_shifted_by_holiday": int(base_exdate != base_third_friday),
        "quote_date_shifted_by_holiday": int(quote_date != raw_quote_date),
    }


# ============================================================
# 3. Generate target table
# ============================================================

def make_target_table(ttm_label, months_ahead):
    rows = []

    if END_BY == "quote_month":

        quote_months = pd.period_range(
            START_QUOTE_MONTH,
            END_QUOTE_MONTH,
            freq="M"
        )

    elif END_BY == "expiration_month":

        start_quote_month = pd.Period(START_QUOTE_MONTH, freq="M")
        end_exp_month = pd.Period(END_EXP_MONTH, freq="M")

        # final expiration month = quote_month + months_ahead
        last_quote_month = end_exp_month - months_ahead

        quote_months = pd.period_range(
            start_quote_month,
            last_quote_month,
            freq="M"
        )

    else:
        raise ValueError("END_BY must be either 'quote_month' or 'expiration_month'.")

    for quote_month in quote_months:

        quote_info = target_quote_date_from_quote_month(quote_month)

        exp_month = quote_month + months_ahead
        exdate, third_fri = monthly_spx_settlement_date(exp_month)

        quote_date = quote_info["quote_date"]

        rows.append({
            "date": quote_date.strftime("%Y%m%d"),
            "exdate": exdate.strftime("%Y%m%d"),

            # audit columns
            "ttm_label": ttm_label,
            "months_ahead": months_ahead,
            "quote_month": quote_month.strftime("%Y%m"),
            "exp_month": exp_month.strftime("%Y%m"),
            "third_friday": third_fri.strftime("%Y%m%d"),
            "raw_quote_date": quote_info["raw_quote_date"].strftime("%Y%m%d"),
            "base_exp_month": quote_info["base_exp_month"].strftime("%Y%m"),
            "base_exdate": quote_info["base_exdate"].strftime("%Y%m%d"),
            "base_third_friday": quote_info["base_third_friday"].strftime("%Y%m%d"),
            "base_exdate_shifted_by_holiday": quote_info["base_exdate_shifted_by_holiday"],
            "exdate_shifted_by_holiday": int(exdate != third_fri),
            "quote_date_shifted_by_holiday": quote_info["quote_date_shifted_by_holiday"],
            "calendar_days_to_exdate": (exdate - quote_date).days,
            "ttm_after_am_settlement": (exdate - quote_date).days - 1,
        })

    target = pd.DataFrame(rows)

    return target


# ============================================================
# 3A. Raw option availability for 90/180-day targets
# ============================================================

def effective_settlement_date(raw_exdate):
    """Convert OptionMetrics raw exdate to the effective settlement date."""
    raw_exdate = to_date(raw_exdate)

    # Historical standard monthly SPX exdates are Saturday-coded.
    if raw_exdate.dayofweek == 5:
        raw_exdate -= pd.Timedelta(days=1)

    # Handles Good Friday and any other exchange holiday.
    return previous_trading_day(raw_exdate)


def option_file_sort_key(path):
    match = re.fullmatch(r"OP(\d{4})_(\d{1,2})\.txt", path.name)
    return int(match.group(1)), int(match.group(2))


def build_option_chain_metrics():
    """
    Build date-by-expiration availability and conservative data sufficiency
    diagnostics from the raw 15-column option files.

    The preliminary screen mirrors the observable part of the MATLAB filters:
    bid > 3/8, ask > bid, finite IV, and the OTM/ATM strike region.  Requiring
    at least 14 observations and 10 unique strikes gives a buffer above the
    MATLAB spline minimum of 10 observations and 6 unique strikes.
    """
    option_files = [
        path
        for path in RAW_OPTION_DIR.glob("OP*.txt")
        if re.fullmatch(r"OP\d{4}_\d{1,2}\.txt", path.name)
    ]
    option_files = sorted(option_files, key=option_file_sort_key)

    if not option_files:
        raise FileNotFoundError(
            f"No monthly option files found in {RAW_OPTION_DIR}"
        )

    all_metrics = []

    for file_number, path in enumerate(option_files, start=1):
        raw = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            usecols=[0, 1, 2, 3, 4, 5, 6, 7, 10],
            names=[
                "secid",
                "date",
                "ttm",
                "cp",
                "k",
                "s",
                "bid",
                "ask",
                "iv",
            ],
        )

        raw = raw.loc[raw["secid"].eq(108105)].copy()
        if raw.empty:
            continue

        raw["date_dt"] = pd.to_datetime(
            raw["date"].astype(str), format="%Y%m%d"
        )
        raw["raw_exdate_dt"] = (
            raw["date_dt"] + pd.to_timedelta(raw["ttm"], unit="D")
        )

        exdate_map = {
            value: effective_settlement_date(value)
            for value in raw["raw_exdate_dt"].drop_duplicates()
        }
        raw["exdate_dt"] = raw["raw_exdate_dt"].map(exdate_map)

        raw["basic_ok"] = (
            raw["bid"].gt(3.0 / 8.0)
            & raw["ask"].gt(raw["bid"])
            & np.isfinite(raw["iv"])
        )
        raw["rough_otm"] = (
            (raw["cp"].eq(1) & raw["k"].ge(raw["s"] - 20))
            | (raw["cp"].eq(2) & raw["k"].le(raw["s"] + 20))
        )
        raw["prelim_ok"] = raw["basic_ok"] & raw["rough_otm"]

        metrics = (
            raw.groupby(["date_dt", "exdate_dt"], as_index=False)
            .agg(
                raw_obs=("k", "size"),
                prelim_obs=("prelim_ok", "sum"),
            )
        )
        unique_k = (
            raw.loc[raw["prelim_ok"]]
            .groupby(["date_dt", "exdate_dt"])["k"]
            .nunique()
            .rename("prelim_unique_k")
            .reset_index()
        )
        metrics = metrics.merge(
            unique_k,
            on=["date_dt", "exdate_dt"],
            how="left",
        )
        metrics["prelim_unique_k"] = (
            metrics["prelim_unique_k"].fillna(0).astype(int)
        )
        all_metrics.append(metrics)

        if file_number % 25 == 0 or file_number == len(option_files):
            print(
                f"Scanned option files: {file_number}/{len(option_files)}",
                flush=True,
            )

    return pd.concat(all_metrics, ignore_index=True)


def select_data_driven_targets(target, ttm_label, chain_metrics):
    """
    Select one sufficiently populated standard-monthly chain per target month.

    Score = absolute maturity error + 2 * absolute quote-date shift.
    Thus an old quote date is retained unless moving a few days materially
    improves maturity or is necessary for contract availability/data quality.
    """
    lower_ttm, upper_ttm = TTM_BOUNDS[ttm_label]
    selected = target.copy()

    audit_rows = []

    for row_number, row in target.iterrows():
        original_date = pd.to_datetime(row["date"], format="%Y%m%d")
        original_exdate = pd.to_datetime(row["exdate"], format="%Y%m%d")

        candidates = chain_metrics.copy()
        candidates["date_shift_days"] = (
            candidates["date_dt"] - original_date
        ).dt.days
        candidates["date_shift_abs"] = candidates["date_shift_days"].abs()
        candidates["final_ttm"] = (
            candidates["exdate_dt"] - candidates["date_dt"]
        ).dt.days - 1

        candidates = candidates.loc[
            candidates["date_shift_abs"].le(MAX_QUOTE_DATE_SHIFT_DAYS)
            & candidates["final_ttm"].between(lower_ttm, upper_ttm)
            & candidates["prelim_obs"].ge(MIN_PRELIMINARY_OBS)
            & candidates["prelim_unique_k"].ge(MIN_PRELIMINARY_UNIQUE_K)
        ].copy()

        if candidates.empty:
            raise RuntimeError(
                "No sufficiently populated option chain found for "
                f"TTM={ttm_label}, original date={row['date']}, "
                f"within +/-{MAX_QUOTE_DATE_SHIFT_DAYS} calendar days."
            )

        candidates["ttm_error_days"] = (
            candidates["final_ttm"] - ttm_label
        ).abs()
        candidates["same_original_pair"] = (
            candidates["date_dt"].eq(original_date)
            & candidates["exdate_dt"].eq(original_exdate)
        )
        candidates["selection_score"] = (
            candidates["ttm_error_days"]
            + QUOTE_DATE_SHIFT_PENALTY * candidates["date_shift_abs"]
        )

        candidates = candidates.sort_values(
            [
                "selection_score",
                "date_shift_abs",
                "ttm_error_days",
                "same_original_pair",
                "prelim_unique_k",
                "prelim_obs",
                "date_dt",
                "exdate_dt",
            ],
            ascending=[True, True, True, False, False, False, True, True],
        )
        choice = candidates.iloc[0]

        chosen_date = choice["date_dt"]
        chosen_exdate = choice["exdate_dt"]
        date_shift_days = int(choice["date_shift_days"])

        if chosen_date == original_date and chosen_exdate == original_exdate:
            selection_mode = "original_pair"
        elif chosen_exdate == original_exdate:
            selection_mode = "nearby_date_same_expiration"
        elif chosen_date == original_date:
            selection_mode = "same_date_nearest_expiration"
        else:
            selection_mode = "nearby_date_nearest_expiration"

        selected.at[row_number, "date"] = chosen_date.strftime("%Y%m%d")
        selected.at[row_number, "exdate"] = chosen_exdate.strftime("%Y%m%d")
        selected.at[row_number, "exp_month"] = chosen_exdate.strftime("%Y%m")

        chosen_exp_period = chosen_exdate.to_period("M")
        chosen_third_friday = third_friday(
            chosen_exp_period.year,
            chosen_exp_period.month,
        )
        selected.at[row_number, "third_friday"] = (
            chosen_third_friday.strftime("%Y%m%d")
        )
        selected.at[row_number, "exdate_shifted_by_holiday"] = int(
            chosen_exdate != chosen_third_friday
        )
        selected.at[row_number, "calendar_days_to_exdate"] = (
            chosen_exdate - chosen_date
        ).days
        selected.at[row_number, "ttm_after_am_settlement"] = int(
            choice["final_ttm"]
        )

        quote_period = pd.Period(row["quote_month"], freq="M")
        selected_months_ahead = (
            (chosen_exp_period.year - quote_period.year) * 12
            + chosen_exp_period.month
            - quote_period.month
        )

        audit_rows.append(
            {
                "original_date": original_date.strftime("%Y%m%d"),
                "original_exdate": original_exdate.strftime("%Y%m%d"),
                "date_shift_days": date_shift_days,
                "selected_months_ahead": selected_months_ahead,
                "ttm_error_days": int(choice["ttm_error_days"]),
                "selection_score": int(choice["selection_score"]),
                "selection_mode": selection_mode,
                "raw_chain_obs": int(choice["raw_obs"]),
                "preliminary_obs": int(choice["prelim_obs"]),
                "preliminary_unique_k": int(choice["prelim_unique_k"]),
                "contract_available": 1,
            }
        )

    audit = pd.DataFrame(audit_rows, index=selected.index)
    selected = pd.concat([selected, audit], axis=1)

    return selected


# ============================================================
# 4. Validation
# ============================================================

def validate_target_table(target, ttm_label):
    date_dt = pd.to_datetime(target["date"], format="%Y%m%d")
    exdate_dt = pd.to_datetime(target["exdate"], format="%Y%m%d")

    duplicated_date = target["date"].duplicated().sum()
    duplicated_pair = target[["date", "exdate"]].duplicated().sum()

    non_trading_dates = [
        d.strftime("%Y%m%d")
        for d in date_dt
        if not is_trading_day(d)
    ]

    non_trading_exdates = [
        d.strftime("%Y%m%d")
        for d in exdate_dt
        if not is_trading_day(d)
    ]

    if duplicated_date > 0:
        print(f"[Warning] TTM_{ttm_label}: duplicated date count = {duplicated_date}")

    if duplicated_pair > 0:
        print(f"[Warning] TTM_{ttm_label}: duplicated date-exdate pair count = {duplicated_pair}")

    if len(non_trading_dates) > 0:
        print(f"[Warning] TTM_{ttm_label}: non-trading quote dates:")
        print(non_trading_dates)

    if len(non_trading_exdates) > 0:
        print(f"[Warning] TTM_{ttm_label}: non-trading exdates:")
        print(non_trading_exdates)

    if len(non_trading_dates) == 0 and len(non_trading_exdates) == 0:
        print(f"[OK] TTM_{ttm_label}: all date/exdate are NYSE trading days.")

    if ttm_label in DATA_DRIVEN_TTM_LABELS:
        lower_ttm, upper_ttm = TTM_BOUNDS[ttm_label]
        final_ttm = (exdate_dt - date_dt).dt.days - 1

        validation_errors = []
        if duplicated_date > 0 or duplicated_pair > 0:
            validation_errors.append("duplicate target dates or pairs")
        if non_trading_dates or non_trading_exdates:
            validation_errors.append("non-trading dates")
        if not final_ttm.between(lower_ttm, upper_ttm).all():
            validation_errors.append("TTM outside configured bounds")
        if target["date_shift_days"].abs().max() > MAX_QUOTE_DATE_SHIFT_DAYS:
            validation_errors.append("quote-date shift outside configured window")
        if not target["contract_available"].eq(1).all():
            validation_errors.append("unavailable option contract")
        if not target["preliminary_obs"].ge(MIN_PRELIMINARY_OBS).all():
            validation_errors.append("insufficient preliminary observations")
        if not target["preliminary_unique_k"].ge(
            MIN_PRELIMINARY_UNIQUE_K
        ).all():
            validation_errors.append("insufficient preliminary unique strikes")

        if validation_errors:
            raise RuntimeError(
                f"TTM_{ttm_label} validation failed: "
                + "; ".join(validation_errors)
            )

        print(
            f"[OK] TTM_{ttm_label}: all {len(target)} months have an available "
            f"chain; final TTM range = {final_ttm.min()}-{final_ttm.max()} days; "
            f"maximum quote-date shift = {target['date_shift_days'].abs().max()} days."
        )


# ============================================================
# 5. Output
# ============================================================

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

chain_metrics = build_option_chain_metrics()
targets_by_label = {}

for ttm_label, months_ahead in TTM_MONTH_MAP.items():
    target = make_target_table(ttm_label, months_ahead)

    if ttm_label in DATA_DRIVEN_TTM_LABELS:
        target = select_data_driven_targets(
            target,
            ttm_label,
            chain_metrics,
        )

    validate_target_table(target, ttm_label)

    targets_by_label[ttm_label] = target

del chain_metrics

# Write only after every target table has passed validation.
all_audit = []

for ttm_label, target in targets_by_label.items():

    simple_path = OUTPUT_DIR / f"TTM_{ttm_label}.csv"
    target[["date", "exdate"]].to_csv(simple_path, index=False)

    print(f"Saved: {simple_path}")
    print(target[["date", "exdate", "calendar_days_to_exdate", "ttm_after_am_settlement"]].tail(5))
    print()

    if SAVE_AUDIT:
        audit_path = OUTPUT_DIR / f"TTM_{ttm_label}_audit.csv"
        target.to_csv(audit_path, index=False)

        print(f"Saved audit: {audit_path}")
        print()

        all_audit.append(target)

if SAVE_AUDIT and len(all_audit) > 0:
    all_audit = pd.concat(all_audit, ignore_index=True)
    all_audit_path = OUTPUT_DIR / "TTM_All_audit.csv"
    all_audit.to_csv(all_audit_path, index=False)

    print(f"Saved combined audit: {all_audit_path}")
