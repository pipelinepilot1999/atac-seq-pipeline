import pandas as pd
from pathlib import Path

metrics_files = snakemake.input.metrics
peaks_files   = snakemake.input.peaks
output_file   = snakemake.output.summary

rows = []

for metrics_file, peaks_file in zip(sorted(metrics_files), sorted(peaks_files)):
    sample = Path(metrics_file).stem.replace(".metrics", "")

    # Parse Picard duplication rate
    dup_rate = "NA"
    with open(metrics_file) as f:
        for line in f:
            if line.startswith("ATAC"):
                fields = line.strip().split("\t")
                if len(fields) >= 9:
                    try:
                        dup_rate = float(fields[8])
                    except:
                        dup_rate = "NA"
                break

    # Count peaks
    peak_count = 0
    with open(peaks_file) as f:
        for line in f:
            if line.strip():
                peak_count += 1

    rows.append({
        "sample": sample,
        "duplication_rate": dup_rate,
        "peak_count": peak_count
    })

df = pd.DataFrame(rows)
df.to_csv(output_file, sep="\t", index=False)
print(df.to_string(index=False))
