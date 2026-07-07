import pandas as pd
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