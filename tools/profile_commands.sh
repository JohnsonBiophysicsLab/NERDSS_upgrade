#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Print NERDSS profiler command blocks without executing them.

Usage:
  tools/profile_commands.sh [options]

Options:
  --mode serial|mpi        Executable mode to profile. Default: serial.
  --profiler TOOL          all, gprof, perf, instruments, gperftools, sample.
                           Default: all.
  --case DIR               Case directory containing the input and .mol files.
                           Default: sample_inputs/VALIDATE_SUITE/homoTrimer.
  --input FILE             Input file name inside --case. Default: parmTri6.inp.
  --seed N                 Fixed seed for repeatable profiling. Default: 12345.
  --ranks N                MPI rank count for --mode mpi. Default: 4.
  --out DIR                Output root for profile runs. Default: profile-runs.
  --help                   Show this help.

Examples:
  tools/profile_commands.sh --mode serial --profiler perf
  tools/profile_commands.sh --mode mpi --profiler gperftools --ranks 4
USAGE
}

mode="serial"
profiler="all"
case_dir="sample_inputs/VALIDATE_SUITE/homoTrimer"
input_file="parmTri6.inp"
seed="12345"
ranks="4"
out_root="profile-runs"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      mode="${2:?missing value for --mode}"
      shift 2
      ;;
    --profiler)
      profiler="${2:?missing value for --profiler}"
      shift 2
      ;;
    --case)
      case_dir="${2:?missing value for --case}"
      shift 2
      ;;
    --input)
      input_file="${2:?missing value for --input}"
      shift 2
      ;;
    --seed)
      seed="${2:?missing value for --seed}"
      shift 2
      ;;
    --ranks)
      ranks="${2:?missing value for --ranks}"
      shift 2
      ;;
    --out)
      out_root="${2:?missing value for --out}"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "$mode" != "serial" && "$mode" != "mpi" ]]; then
  echo "--mode must be serial or mpi" >&2
  exit 2
fi

case "$profiler" in
  all|gprof|perf|instruments|gperftools|sample) ;;
  *)
    echo "--profiler must be all, gprof, perf, instruments, gperftools, or sample" >&2
    exit 2
    ;;
esac

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
case_name="$(basename "$case_dir")"
run_root="$repo_root/$out_root/${mode}-${case_name}"
exe_name="nerdss"
make_target="serial"
runner=""
make_override=""

if [[ "$mode" == "mpi" ]]; then
  exe_name="nerdss_mpi"
  make_target="mpi"
  runner="mpirun -np $ranks "
  make_override=" CC=mpicxx"
fi

print_header() {
  local title="$1"
  printf '\n# %s\n' "$title"
}

print_prep() {
  local run_dir="$1"
  printf 'cd %q\n' "$repo_root"
  printf 'make clean\n'
  printf 'make %s%s%s\n' "$make_target" "$2" "$make_override"
  printf 'mkdir -p %q\n' "$run_dir"
  printf 'cp %q/*.inp %q/*.mol %q/ 2>/dev/null || true\n' "$repo_root/$case_dir" "$repo_root/$case_dir" "$run_dir"
}

print_run_command() {
  local run_dir="$1"
  local prefix="${2:-}"
  printf '( cd %q && %s%s%q -f %q -s %q )\n' \
    "$run_dir" "$prefix" "$runner" "$repo_root/bin/$exe_name" "$input_file" "$seed"
}

want() {
  [[ "$profiler" == "all" || "$profiler" == "$1" ]]
}

cat <<EOF
# NERDSS profile command guide
# Repository: $repo_root
# Mode: $mode
# Case: $case_dir
# Input: $input_file
# Seed: $seed
EOF

if want gprof; then
  run_dir="$run_root-gprof"
  print_header "gprof"
  print_prep "$run_dir" " profile"
  if [[ "$mode" == "mpi" ]]; then
    print_run_command "$run_dir" "GMON_OUT_PREFIX=$run_dir/gmon "
    printf 'for f in %q/gmon.*; do gprof %q "$f" > "$f.txt"; done\n' "$run_dir" "$repo_root/bin/$exe_name"
  else
    print_run_command "$run_dir"
    printf 'gprof %q %q/gmon.out > %q/gprof.txt\n' "$repo_root/bin/$exe_name" "$run_dir" "$run_dir"
  fi
fi

if want perf; then
  run_dir="$run_root-perf"
  print_header "Linux perf"
  print_prep "$run_dir" ""
  if [[ "$mode" == "mpi" ]]; then
    printf '( cd %q && perf record -F 99 -g --call-graph dwarf -o perf.data -- mpirun -np %q %q -f %q -s %q )\n' \
      "$run_dir" "$ranks" "$repo_root/bin/$exe_name" "$input_file" "$seed"
  else
    printf '( cd %q && perf record -F 99 -g --call-graph dwarf -o perf.data -- %q -f %q -s %q )\n' \
      "$run_dir" "$repo_root/bin/$exe_name" "$input_file" "$seed"
  fi
  printf 'perf report -i %q/perf.data\n' "$run_dir"
  printf 'perf script -i %q/perf.data > %q/perf.script.txt\n' "$run_dir" "$run_dir"
fi

if want instruments; then
  run_dir="$run_root-instruments"
  print_header "macOS Instruments"
  print_prep "$run_dir" ""
  if [[ "$mode" == "mpi" ]]; then
    printf '( cd %q && xcrun xctrace record --template "Time Profiler" --output TimeProfiler.trace --launch -- mpirun -np %q %q -f %q -s %q )\n' \
      "$run_dir" "$ranks" "$repo_root/bin/$exe_name" "$input_file" "$seed"
  else
    printf '( cd %q && xcrun xctrace record --template "Time Profiler" --output TimeProfiler.trace --launch -- %q -f %q -s %q )\n' \
      "$run_dir" "$repo_root/bin/$exe_name" "$input_file" "$seed"
  fi
  printf 'xcrun xctrace export --input %q/TimeProfiler.trace --xpath '"'"'/trace-toc/run/data/table'"'"' > %q/xctrace-table.xml\n' "$run_dir" "$run_dir"
  printf 'open %q/TimeProfiler.trace\n' "$run_dir"
fi

if want gperftools; then
  run_dir="$run_root-gperftools"
  print_header "gperftools"
  print_prep "$run_dir" " profile"
  if [[ "$mode" == "mpi" ]]; then
    print_run_command "$run_dir"
    printf 'for f in %q/profile_output_*.prof; do pprof --text %q "$f" > "${f%%.prof}.txt"; done\n' "$run_dir" "$repo_root/bin/$exe_name"
  else
    print_run_command "$run_dir" "CPUPROFILE=$run_dir/cpu.prof "
    printf 'pprof --text %q %q/cpu.prof > %q/pprof.txt\n' "$repo_root/bin/$exe_name" "$run_dir" "$run_dir"
  fi
fi

if want sample; then
  print_header "macOS sample"
  cat <<'EOF'
# Run NERDSS in one terminal, then sample it from another:
pgrep -fl nerdss
sample <PID> 10 -file sample.txt
EOF
fi
