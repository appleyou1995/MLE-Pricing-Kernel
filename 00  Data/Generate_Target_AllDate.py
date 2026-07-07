import pandas as pd
from pathlib import Path
from calendar import monthcalendar, FRIDAY
import exchange_calendars as xcals


# ============================================================
# 0. Parameters
# ============================================================

START_QUOTE_MONTH = "1996-01"
END_QUOTE_MONTH   = "2025-08"

END_BY = "quote_month"
END_EXP_MONTH = "2025-08"   # 這行在 END_BY = "quote_month" 時不會用到

OUTPUT_DIR = Path(
    r"D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel\Code\01  Output"
)

SAVE_AUDIT = True

# TTM label 對應到幾個月後的 traditional monthly SPX expiration
# 注意：這裡的 60/90/180 是沿用舊檔命名，不是嚴格 calendar days。
TTM_MONTH_MAP = {
    30: 1,
    60: 2,
    90: 3,
    180: 6,
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


# ============================================================
# 5. Output
# ============================================================

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

all_audit = []

for ttm_label, months_ahead in TTM_MONTH_MAP.items():

    target = make_target_table(ttm_label, months_ahead)

    validate_target_table(target, ttm_label)

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