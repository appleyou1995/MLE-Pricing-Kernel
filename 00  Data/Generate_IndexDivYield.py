import pandas as pd
from pathlib import Path


# ============================================================
# 0. 路徑設定
# ============================================================

data_dir = Path(
    r"D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel\Data"
)

output_dir = Path(
    r"D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel\Code\01  Output"
)

input_file = data_dir / "raw_IndexDivYield_19962025.csv"
output_file = output_dir / "IndexDivYield19962025.txt"


# ============================================================
# 1. 讀取 csv
# ============================================================

df = pd.read_csv(input_file)

# 欄位名稱保險處理：去掉前後空白
df.columns = df.columns.str.strip()


# ============================================================
# 2. 篩選 secid = 108105
# ============================================================

df = df[df["secid"] == 108105].copy()


# ============================================================
# 3. 日期由小排到大
# ============================================================

df["date"] = pd.to_numeric(df["date"], errors="coerce")
df["rate"] = pd.to_numeric(df["rate"], errors="coerce")

df = df.dropna(subset=["date", "rate"])
df["date"] = df["date"].astype(int)

df = df.sort_values("date")


# ============================================================
# 4. rate 除以 100
# ============================================================

df["rate"] = df["rate"] / 100


# ============================================================
# 5. 輸出成 txt：無標題列、無 index、空白分隔
# ============================================================

output_dir.mkdir(parents=True, exist_ok=True)

with open(output_file, "w", encoding="utf-8", newline="\n") as f:
    for _, row in df.iterrows():
        secid = int(row["secid"])
        date = int(row["date"])
        rate = format(row["rate"], ".10g")   # 避免多餘尾數 0，格式接近舊檔

        f.write(f"{secid} {date} {rate}\n")


# ============================================================
# 6. 檢查輸出
# ============================================================

print(f"Saved: {output_file}")
print(f"Number of rows: {len(df):,}")
print(df.head())
print(df.tail())