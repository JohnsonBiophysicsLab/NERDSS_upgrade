#!/usr/bin/env bash
#
# run_profile.sh — actually run a per-function timing profile of NERDSS.
#
# Unlike tools/profile_commands.sh (which only prints command blocks), this
# script sets up an isolated run directory, applies an optional nItr cap, runs
# the chosen profiler, and prints a per-function breakdown.
#
#   --profiler instruments  Instruments Time Profiler (xctrace) -> per-function
#                           self/inclusive SECONDS via tools/profile_aggregate.py (default)
#   --profiler pprof        macOS call graph: /usr/bin/sample -> instrumentsToPprof
#                           -> `go tool pprof` (real symbolized graph; needs a long
#                           enough run, so raise --nitr / --sample-seconds)
#   --profiler perf         Linux: perf record -g --call-graph dwarf -> per-function
#                           seconds; converts perf.data -> pprof (if perf_to_profile
#                           present) for the same -peek/-png/-http workflow
#   --profiler gperftools   libprofiler CPU profile via preload (Linux; unwinding is
#                           broken on arm64 macOS -- use pprof mode there instead)
#   --profiler both         macOS: instruments + pprof ; Linux: perf
#
# Examples:
#   tools/run_profile.sh --case sample_inputs/VALIDATE_SUITE/homoTrimer \
#       --input parmTri6.inp --nitr 2000
#
#   # macOS call graph:
#   tools/run_profile.sh --case sample_inputs/VALIDATE_SUITE/clathrin \
#       --input parms_clath_kon1uM.inp --nitr 300000 --profiler pprof --build
#
#   # Linux perf:
#   tools/run_profile.sh --case sample_inputs/VALIDATE_SUITE/clathrin \
#       --input parms_clath_kon1uM.inp --nitr 100000 --profiler perf --build
#
set -euo pipefail

# ---------------------------------------------------------------- defaults
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

case_dir=""
input_file=""
seed="12345"
nitr=""                       # empty => run input unchanged
profiler="instruments"        # instruments | gperftools | both
bin="$repo_root/bin/nerdss"
out_root="$repo_root/profile-runs"
top="25"
do_build="0"
freq="1000"                   # gperftools sampling Hz
sample_seconds="15"           # pprof mode: how long to /usr/bin/sample

usage() { awk 'NR>1 && /^#/{sub(/^# ?/,""); print; next} NR>1{exit}' "${BASH_SOURCE[0]}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case)      case_dir="$2"; shift 2 ;;
    --input)     input_file="$2"; shift 2 ;;
    --seed)      seed="$2"; shift 2 ;;
    --nitr)      nitr="$2"; shift 2 ;;
    --profiler)  profiler="$2"; shift 2 ;;
    --bin)       bin="$2"; shift 2 ;;
    --out)       out_root="$2"; shift 2 ;;
    --top)       top="$2"; shift 2 ;;
    --freq)      freq="$2"; shift 2 ;;
    --sample-seconds) sample_seconds="$2"; shift 2 ;;
    --build)     do_build="1"; shift ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$case_dir"   ]] || { echo "ERROR: --case is required"  >&2; exit 2; }
[[ -n "$input_file" ]] || { echo "ERROR: --input is required" >&2; exit 2; }
case "$profiler" in instruments|pprof|perf|gperftools|both) ;; *)
  echo "ERROR: --profiler must be instruments, pprof, perf, gperftools, or both" >&2; exit 2 ;;
esac

# case_dir may be relative to repo root
[[ "$case_dir" = /* ]] || case_dir="$repo_root/$case_dir"
[[ -d "$case_dir" ]] || { echo "ERROR: case dir not found: $case_dir" >&2; exit 2; }

# ---------------------------------------------------------------- build
if [[ "$do_build" == "1" ]]; then
  echo ">> Building optimized serial binary (make serial)"
  ( cd "$repo_root" && make serial )
fi
[[ -x "$bin" ]] || { echo "ERROR: binary not found/executable: $bin (use --build)" >&2; exit 2; }

# ---------------------------------------------------------------- run dir
case_name="$(basename "$case_dir")"
run_dir="$out_root/${profiler}-${case_name}"
mkdir -p "$run_dir"
echo ">> Run dir: $run_dir"
cp "$case_dir"/*.inp "$run_dir"/ 2>/dev/null || true
cp "$case_dir"/*.mol "$run_dir"/ 2>/dev/null || true
[[ -f "$run_dir/$input_file" ]] || { echo "ERROR: input not in case dir: $input_file" >&2; exit 2; }

# Apply nItr cap in the isolated copy only (raw inputs often run ~forever).
if [[ -n "$nitr" ]]; then
  # matches "nItr = <val>" case-insensitively, preserves leading whitespace
  perl -0pi -e "s/^(\s*nItr\s*=\s*).*$/\${1}$nitr/mi" "$run_dir/$input_file"
  echo ">> Capped nItr = $nitr in $run_dir/$input_file"
fi

run_cmd=( "$bin" -f "$input_file" -s "$seed" )
echo ">> Command: ${run_cmd[*]} (seed $seed)"

# Measure real CPU seconds (user+sys) from a plain, un-profiled run so the
# aggregator can scale sample shares to real time. Uses `/usr/bin/time -p`
# (POSIX; on Linux this is the `time` package). Returns "" if unavailable.
measure_cpu_seconds() {
  [[ -x /usr/bin/time ]] || { echo ""; return; }
  local tlog; tlog="$(mktemp)"
  ( cd "$run_dir" && /usr/bin/time -p "${run_cmd[@]}" >/dev/null 2>"$tlog" ) || true
  rm -rf "$run_dir/DATA" "$run_dir/RESTARTS"
  awk '/^user/{u=$2} /^sys/{s=$2} END{if(u!=""||s!="")printf "%.2f", u+s}' "$tlog"
  rm -f "$tlog"
}
echo ">> Measuring CPU time (plain run) ..."
cpu_seconds="$(measure_cpu_seconds)"
if [[ -n "$cpu_seconds" ]]; then
  echo ">> Measured CPU time: ${cpu_seconds}s (user+sys)"
else
  echo ">> /usr/bin/time not available; seconds columns will be sampling-derived"
fi
# aggregator flag: only scale to real seconds when we measured them
secs_flag=()
[[ -n "$cpu_seconds" ]] && secs_flag=(--total-seconds "$cpu_seconds")

# ---------------------------------------------------------------- instruments
run_instruments() {
  command -v xcrun >/dev/null || { echo "ERROR: xcrun not found (macOS only)" >&2; return 2; }
  local trace="$run_dir/TimeProfiler.trace"
  local xml="$run_dir/time-profile.xml"
  rm -rf "$trace" "$run_dir/DATA" "$run_dir/RESTARTS"
  echo ">> [instruments] recording Time Profiler ..."
  ( cd "$run_dir" && xcrun xctrace record --template "Time Profiler" \
      --output "$trace" --launch -- "${run_cmd[@]}" >/dev/null )
  echo ">> [instruments] exporting time-profile table ..."
  xcrun xctrace export --input "$trace" \
    --xpath '/trace-toc/run[@number="1"]/data/table[@schema="time-profile"]' > "$xml"
  echo "----------------------------------------------------------------------"
  python3 "$script_dir/profile_aggregate.py" "$xml" --top "$top" "${secs_flag[@]}"
  echo "----------------------------------------------------------------------"
  echo ">> trace: $trace"
  echo ">> xml:   $xml"
}

# ---------------------------------------------------------------- pprof (macOS)
# Real symbolized pprof call graph on macOS: /usr/bin/sample the running process,
# then convert with google/instrumentsToPprof (--format=sample). This works where
# gperftools' arm64 unwinder does not. Needs a run long enough to sample, so pass
# a larger --nitr (and/or --sample-seconds).
run_pprof() {
  command -v go >/dev/null || { echo "ERROR: 'go' not found (go tool pprof)" >&2; return 2; }
  local itp; itp="$(command -v instrumentsToPprof 2>/dev/null || echo "$(go env GOPATH)/bin/instrumentsToPprof")"
  [[ -x "$itp" ]] || { echo "ERROR: instrumentsToPprof not found. Install: go install github.com/google/instrumentsToPprof@latest" >&2; return 2; }
  [[ -x /usr/bin/sample ]] || { echo "ERROR: /usr/bin/sample not found (macOS only)" >&2; return 2; }
  local sfile="$run_dir/sample.txt" pb="$run_dir/cpu.pb.gz"
  rm -f "$sfile" "$pb"; rm -rf "$run_dir/DATA" "$run_dir/RESTARTS"
  echo ">> [pprof] launching + sampling for ${sample_seconds}s ..."
  pushd "$run_dir" >/dev/null
  "${run_cmd[@]}" >/dev/null 2>&1 &
  local pid=$!
  /usr/bin/sample "$pid" "$sample_seconds" -file "$sfile" >/dev/null 2>&1 || true
  wait "$pid" 2>/dev/null || true
  popd >/dev/null
  rm -rf "$run_dir/DATA" "$run_dir/RESTARTS"
  [[ -s "$sfile" ]] || { echo "ERROR: sample produced no data (did the run exit too fast? raise --nitr)" >&2; return 2; }
  "$itp" --format=sample -output "$pb" "$sfile"
  echo "----------------------------------------------------------------------"
  echo "## pprof flat (top $top)"
  go tool pprof -nodecount="$top" -top "$pb" 2>/dev/null || true
  echo
  echo "## pprof call tree (top $top)"
  go tool pprof -nodecount="$top" -tree "$pb" 2>/dev/null || true
  echo "----------------------------------------------------------------------"
  echo ">> pprof profile: $pb"
  echo ">> callers of X:  go tool pprof -peek=X $pb"
  echo ">> web graph:     go tool pprof -http=: $pb   (needs graphviz: brew install graphviz)"
}

# ---------------------------------------------------------------- gperftools
find_libprofiler() {
  local names=(libprofiler.dylib libprofiler.so)
  local dirs=()
  if command -v brew >/dev/null; then
    local cellar; cellar="$(brew --cellar gperftools 2>/dev/null || true)"
    [[ -n "$cellar" ]] && dirs+=( "$cellar"/*/lib )
    dirs+=( "$(brew --prefix 2>/dev/null)/lib" )
  fi
  dirs+=( /usr/local/lib /usr/lib /usr/lib/x86_64-linux-gnu /opt/homebrew/lib )
  for d in "${dirs[@]}"; do
    for n in "${names[@]}"; do
      [[ -f "$d/$n" ]] && { echo "$d/$n"; return 0; }
    done
  done
  return 1
}

run_gperftools() {
  local lp; lp="$(find_libprofiler || true)"
  [[ -n "$lp" ]] || { echo "ERROR: libprofiler not found (brew install gperftools / apt install libgoogle-perftools-dev)" >&2; return 2; }
  command -v go >/dev/null || { echo "ERROR: 'go' not found for 'go tool pprof'. Install go or gperftools' pprof." >&2; return 2; }
  local prof="$run_dir/cpu.prof"
  local preload="LD_PRELOAD"
  [[ "$(uname)" == "Darwin" ]] && preload="DYLD_INSERT_LIBRARIES"
  rm -f "$prof"; rm -rf "$run_dir/DATA" "$run_dir/RESTARTS"
  echo ">> [gperftools] using $lp"
  echo ">> [gperftools] recording CPU profile (${freq} Hz) ..."
  ( cd "$run_dir" && env "$preload=$lp" CPUPROFILE="$prof" CPUPROFILE_FREQUENCY="$freq" "${run_cmd[@]}" >/dev/null )
  [[ -s "$prof" ]] || { echo "ERROR: no profile produced ($prof)" >&2; return 2; }
  echo "----------------------------------------------------------------------"
  echo "## pprof flat (top $top)"
  go tool pprof -nodecount="$top" -top "$bin" "$prof" 2>/dev/null || true
  echo
  echo "## pprof call tree (top $top)"
  go tool pprof -nodecount="$top" -tree "$bin" "$prof" 2>/dev/null || true
  echo "----------------------------------------------------------------------"
  echo ">> profile: $prof"
  echo ">> interactive:  go tool pprof $bin $prof"
  echo ">> caller/callee: go tool pprof -peek=<function> $bin $prof"
  echo ">> NOTE: on arm64 macOS gperftools stack unwinding is often broken;"
  echo ">>       prefer Instruments there and use gperftools/pprof on Linux."
}

# ---------------------------------------------------------------- perf (Linux)
# Sample the optimized binary with perf, print per-function overhead (scaled to
# real seconds when CPU time was measured), and — when the converter is present —
# emit a pprof profile so the SAME graph/-peek/-http/-png workflow applies.
run_perf() {
  command -v perf >/dev/null || { echo "ERROR: perf not found (Linux; install linux-tools / linux-perf)" >&2; return 2; }
  local data="$run_dir/perf.data" pb="$run_dir/cpu.pb.gz"
  rm -f "$data" "$pb"; rm -rf "$run_dir/DATA" "$run_dir/RESTARTS"
  echo ">> [perf] recording (-F ${freq} -g --call-graph dwarf) ..."
  ( cd "$run_dir" && perf record -F "$freq" -g --call-graph dwarf -o "$data" -- "${run_cmd[@]}" >/dev/null 2>&1 ) || {
    echo "NOTE: perf record failed; you may need 'sysctl kernel.perf_event_paranoid=-1' or sudo." >&2; return 2; }
  rm -rf "$run_dir/DATA" "$run_dir/RESTARTS"

  echo "----------------------------------------------------------------------"
  echo "## perf per-function (flat, top $top)$( [[ -n "$cpu_seconds" ]] && echo " — sec scaled to ${cpu_seconds}s CPU" )"
  # Flat symbol overhead; append real seconds = overhead% * measured CPU seconds.
  perf report -i "$data" --stdio -g none --percent-limit 0.3 2>/dev/null \
    | awk -v cpu="$cpu_seconds" -v n="$top" '
        /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
        {
          pct=$1; sub(/%$/,"",pct); pct+=0
          sym=$0; sub(/.*\][[:space:]]+/, "", sym)   # strip up to "[.] "/"[k] "
          if (cpu=="") printf "%6.2f%%  %s\n", pct, sym
          else         printf "%8.3f s  %6.2f%%  %s\n", cpu*pct/100, pct, sym
          if (++c>=n) exit
        }'

  echo
  # Unified pprof path: go tool pprof reads perf.data when google's
  # perf_to_profile (perf_data_converter) is installed.
  if command -v perf_to_profile >/dev/null && command -v go >/dev/null; then
    perf_to_profile -i "$data" -o "$pb" 2>/dev/null || true
  fi
  if [[ -s "$pb" ]] && command -v go >/dev/null; then
    echo "## pprof call tree (top $top)"
    go tool pprof -nodecount="$top" -tree "$pb" 2>/dev/null || true
    echo ">> pprof profile: $pb  (same workflow as macOS pprof mode)"
    echo ">> callers:   go tool pprof -peek=<fn> $pb"
    echo ">> graph img: go tool pprof -png -nodefraction=0.01 -output cg.png $pb   (needs graphviz)"
    echo ">> web view:  go tool pprof -http=: $pb"
  else
    echo ">> For a call graph, either:"
    echo ">>   perf report -i $data          # interactive TUI call graph"
    echo ">>   perf script -i $data | stackcollapse-perf.pl | flamegraph.pl > flame.svg   # FlameGraph"
    echo ">>   install google/perf_data_converter (perf_to_profile) to reuse the pprof/-png/-http workflow"
  fi
  echo "----------------------------------------------------------------------"
}

# ---------------------------------------------------------------- dispatch
case "$profiler" in
  instruments) run_instruments ;;
  pprof)       run_pprof ;;
  perf)        run_perf ;;
  gperftools)  run_gperftools ;;
  both)        if [[ "$(uname)" == "Darwin" ]]; then run_instruments; echo; run_pprof;
               else run_perf; fi ;;
esac
