#!/usr/bin/env bash
# Download the AHRQ MEPS public use files used in this study (Stata format, ~70 MB zipped).
set -euo pipefail
cd "$(dirname "$0")/.." && mkdir -p data/raw/meps && cd data/raw/meps
for f in h244 h231 h229g h229f h229a h239a; do
  curl -sSL -A "Mozilla/5.0" -o ${f}dta.zip "https://meps.ahrq.gov/mepsweb/data_files/pufs/$f/${f}dta.zip"
  unzip -oq ${f}dta.zip
done
ls -la *.dta
