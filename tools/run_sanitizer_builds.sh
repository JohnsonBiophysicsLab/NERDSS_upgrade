#!/usr/bin/env bash
set -euo pipefail

build_root="${NERDSS_SANITIZER_BUILD_ROOT:-build/sanitizers}"
jobs="${NERDSS_BUILD_JOBS:-}"

if [ "$#" -eq 0 ]; then
  set -- address undefined combined
fi

run_build() {
  mode="$1"
  build_dir="$build_root/$mode"

  asan=OFF
  ubsan=OFF
  case "$mode" in
    address|asan)
      asan=ON
      ;;
    undefined|ubsan)
      ubsan=ON
      ;;
    combined|both|sanitizers)
      asan=ON
      ubsan=ON
      ;;
    *)
      echo "Unknown sanitizer mode: $mode" >&2
      echo "Expected one of: address, undefined, combined" >&2
      return 2
      ;;
  esac

  cmake \
    -S . \
    -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Debug \
    -DNERDSS_ENABLE_ASAN="$asan" \
    -DNERDSS_ENABLE_UBSAN="$ubsan"

  if [ -n "$jobs" ]; then
    cmake --build "$build_dir" --target nerdss --parallel "$jobs"
  else
    cmake --build "$build_dir" --target nerdss
  fi
}

for mode in "$@"; do
  run_build "$mode"
done
