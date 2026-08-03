"""Build RND target-date CSV files from Hsieh's contract workbook.

The workbook's expiration date is preserved as ``raw_exdate`` in the audit
files.  The two-column RND input files use an effective settlement date:

1. Saturday-coded expirations are moved to Friday.
2. If that date is not an NYSE session, it is moved to the previous session.

This keeps the advisor-provided contract selection while remaining compatible
with Procedure__F2010_RND_MultipleDate.m and the realized-return program.
"""

import os
from pathlib import Path

import pandas as pd


OUTPUT_DIR = Path(__file__).resolve().parent
CONTRACT_FILE = OUTPUT_DIR / "contract_list_2026.xlsx"
PROJECT_ROOT = OUTPUT_DIR.parent.parent
DATA_DIR = Path(os.environ.get("MLE_PROJECT_DATA_DIR", PROJECT_ROOT / "Data"))
SPX_FILE = DATA_DIR / "raw_SPX_19962025.csv"

SHEET_BY_TTM = {
    30: "30D",
    60: "60D",
    90: "90D",
    180: "180",
}

def load_spx_trading_dates() -> pd.DatetimeIndex:
    spx = pd.read_csv(SPX_FILE, usecols=["secid", "date"])
    spx.columns = spx.columns.str.strip().str.lower()
    spx = spx.loc[pd.to_numeric(spx["secid"], errors="coerce").eq(108105)].copy()
    dates = pd.to_datetime(
        pd.to_numeric(spx["date"], errors="coerce").dropna().astype(int).astype(str),
        format="%Y%m%d",
        errors="raise",
    )
    return pd.DatetimeIndex(dates.drop_duplicates().sort_values())


TRADING_DATES = load_spx_trading_dates()
TRADING_DATE_SET = set(TRADING_DATES)


def previous_spx_trading_day(value: pd.Timestamp) -> pd.Timestamp:
    value = pd.Timestamp(value).normalize()
    position = TRADING_DATES.searchsorted(value, side="right") - 1
    if position < 0:
        raise ValueError(f"No SPX trading date available on or before {value.date()}.")
    return pd.Timestamp(TRADING_DATES[position]).normalize()


def effective_exdate(raw_exdate: pd.Timestamp) -> tuple[pd.Timestamp, bool, bool]:
    raw_exdate = pd.Timestamp(raw_exdate).normalize()
    candidate = raw_exdate
    saturday_adjusted = candidate.dayofweek == 5
    if saturday_adjusted:
        candidate -= pd.Timedelta(days=1)

    effective = previous_spx_trading_day(candidate)
    holiday_adjusted = effective != candidate
    return effective, saturday_adjusted, holiday_adjusted


def load_contract_sheet(ttm: int, sheet_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = pd.read_excel(
        CONTRACT_FILE,
        sheet_name=sheet_name,
        usecols=["date", "exdate"],
        engine="openpyxl",
    )
    source.columns = source.columns.str.strip().str.lower()
    source = source.dropna(subset=["date", "exdate"]).copy()
    source["date"] = pd.to_datetime(source["date"], errors="raise").dt.normalize()
    source["raw_exdate"] = pd.to_datetime(source["exdate"], errors="raise").dt.normalize()
    source = source.drop(columns="exdate")
    source = source.sort_values("date").reset_index(drop=True)

    if source["date"].duplicated().any():
        duplicates = source.loc[source["date"].duplicated(False), "date"]
        raise ValueError(f"TTM {ttm} contains duplicated quote dates: {duplicates.tolist()}")

    invalid_quote_dates = [
        value.strftime("%Y-%m-%d")
        for value in source["date"]
        if value not in TRADING_DATE_SET
    ]
    if invalid_quote_dates:
        raise ValueError(f"TTM {ttm} contains non-session quote dates: {invalid_quote_dates}")

    adjustments = source["raw_exdate"].map(effective_exdate)
    source["effective_exdate"] = adjustments.map(lambda value: value[0])
    source["saturday_adjusted"] = adjustments.map(lambda value: value[1])
    source["holiday_adjusted"] = adjustments.map(lambda value: value[2])

    if (source["effective_exdate"] <= source["date"]).any():
        raise ValueError(f"TTM {ttm} contains non-positive effective maturities.")

    main = pd.DataFrame({
        "date": source["date"].dt.strftime("%Y%m%d"),
        "exdate": source["effective_exdate"].dt.strftime("%Y%m%d"),
    })

    audit = pd.DataFrame({
        "ttm_label": ttm,
        "source_sheet": sheet_name,
        "date": source["date"].dt.strftime("%Y%m%d"),
        "raw_exdate": source["raw_exdate"].dt.strftime("%Y%m%d"),
        "exdate": source["effective_exdate"].dt.strftime("%Y%m%d"),
        "raw_calendar_ttm": (source["raw_exdate"] - source["date"]).dt.days,
        "effective_calendar_ttm": (
            source["effective_exdate"] - source["date"]
        ).dt.days,
        "final_ttm_after_am": (
            source["effective_exdate"] - source["date"]
        ).dt.days - 1,
        "raw_exdate_weekday": source["raw_exdate"].dt.day_name(),
        "saturday_adjusted": source["saturday_adjusted"].astype(int),
        "holiday_adjusted": source["holiday_adjusted"].astype(int),
    })
    return main, audit


def main() -> None:
    if not CONTRACT_FILE.exists():
        raise FileNotFoundError(f"Contract workbook not found: {CONTRACT_FILE}")

    all_audits = []
    for ttm, sheet_name in SHEET_BY_TTM.items():
        target, audit = load_contract_sheet(ttm, sheet_name)
        target_file = OUTPUT_DIR / f"TTM_{ttm}.csv"
        audit_file = OUTPUT_DIR / f"TTM_{ttm}_audit.csv"
        target.to_csv(target_file, index=False)
        audit.to_csv(audit_file, index=False)
        all_audits.append(audit)

        print(
            f"TTM {ttm}: rows={len(target)}, "
            f"Saturday adjustments={audit['saturday_adjusted'].sum()}, "
            f"holiday adjustments={audit['holiday_adjusted'].sum()}, "
            f"effective TTM={audit['effective_calendar_ttm'].min()}-"
            f"{audit['effective_calendar_ttm'].max()}"
        )

    combined = pd.concat(all_audits, ignore_index=True)
    combined.to_csv(OUTPUT_DIR / "TTM_All_audit.csv", index=False)
    print(f"Saved target files to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
