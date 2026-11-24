#!/usr/bin/env bash
# Usage: bash .claude/scripts/check-chapter.sh 04-models.Rmd
# Knits one chapter standalone, reports RENDER: PASS/FAIL and unresolved citations.
cd "$(dirname "$0")/../.." || exit 2
ch="$1"
[ -f "$ch" ] || { echo "RENDER: FAIL - no such file $ch"; exit 2; }
export R_LIBS=C:/Users/miken/AppData/Local/R/win-library/4.4
export R_LIBS_USER=$R_LIBS
export RSTUDIO_PANDOC="C:/Program Files/Quarto/bin/tools"
export META_BOOK_LOCAL_CHECK=1
R="/c/Program Files/R/R-4.4.3/bin/Rscript.exe"
log="${TEMP:-/tmp}/render-$(basename "$ch" .Rmd).log"
"$R" .claude/scripts/check-chapter.R "$ch" > "$log" 2>&1
code=$?
bad=$(grep -o 'Citeproc: citation [^ ]* not found' "$log" | sort -u)
if [ -n "$bad" ]; then echo "RENDER: FAIL - unresolved citations:"; echo "$bad"; code=1
elif [ $code -eq 0 ] && grep -q 'Output created' "$log"; then echo "RENDER: PASS ($ch)"
else echo "RENDER: FAIL (exit $code) log: $log"; tail -40 "$log"; fi
exit $code
