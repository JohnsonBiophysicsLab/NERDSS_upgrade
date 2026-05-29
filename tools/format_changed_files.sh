#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: tools/format_changed_files.sh [--check] [--base <ref>] [--all]

Formats C/C++ files touched by the current branch. By default, files are
compared against origin/master when available, otherwise master.

Options:
  --check       fail if any changed C/C++ file is not clang-format clean
  --base <ref> compare against a specific git ref
  --all         format/check every tracked C/C++ source file in rollout scope
  -h, --help    show this help text
USAGE
}

mode="format"
base_ref=""
all_files=0

while (($#)); do
  case "$1" in
    --check)
      mode="check"
      shift
      ;;
    --base)
      if (($# < 2)); then
        echo "error: --base requires a git ref" >&2
        exit 2
      fi
      base_ref="$2"
      shift 2
      ;;
    --all)
      all_files=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "error: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

clang_format_bin="${CLANG_FORMAT:-}"
if [[ -z "$clang_format_bin" ]]; then
  if command -v clang-format >/dev/null 2>&1; then
    clang_format_bin="clang-format"
  elif [[ -x /opt/homebrew/opt/llvm/bin/clang-format ]]; then
    clang_format_bin="/opt/homebrew/opt/llvm/bin/clang-format"
  elif [[ -x /usr/local/opt/llvm/bin/clang-format ]]; then
    clang_format_bin="/usr/local/opt/llvm/bin/clang-format"
  fi
fi

if [[ -z "$clang_format_bin" || ! -x "$(command -v "$clang_format_bin" 2>/dev/null || printf '%s' "$clang_format_bin")" ]]; then
  echo "error: clang-format is required but was not found." >&2
  echo "Set CLANG_FORMAT=/path/to/clang-format or add clang-format to PATH." >&2
  exit 127
fi

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

if [[ -z "$base_ref" ]]; then
  if git rev-parse --verify --quiet origin/master >/dev/null; then
    base_ref="origin/master"
  else
    base_ref="master"
  fi
fi

is_source_file() {
  case "$1" in
    EXEs/*.c|EXEs/*.cc|EXEs/*.cpp|EXEs/*.cxx|EXEs/*.h|EXEs/*.hh|EXEs/*.hpp|EXEs/*.hxx) return 0 ;;
    include/*.c|include/*.cc|include/*.cpp|include/*.cxx|include/*.h|include/*.hh|include/*.hpp|include/*.hxx) return 0 ;;
    include/**/*.c|include/**/*.cc|include/**/*.cpp|include/**/*.cxx|include/**/*.h|include/**/*.hh|include/**/*.hpp|include/**/*.hxx) return 0 ;;
    mpi_proto/*.c|mpi_proto/*.cc|mpi_proto/*.cpp|mpi_proto/*.cxx|mpi_proto/*.h|mpi_proto/*.hh|mpi_proto/*.hpp|mpi_proto/*.hxx) return 0 ;;
    mpi_proto/**/*.c|mpi_proto/**/*.cc|mpi_proto/**/*.cpp|mpi_proto/**/*.cxx|mpi_proto/**/*.h|mpi_proto/**/*.hh|mpi_proto/**/*.hpp|mpi_proto/**/*.hxx) return 0 ;;
    src/*.c|src/*.cc|src/*.cpp|src/*.cxx|src/*.h|src/*.hh|src/*.hpp|src/*.hxx) return 0 ;;
    src/**/*.c|src/**/*.cc|src/**/*.cpp|src/**/*.cxx|src/**/*.h|src/**/*.hh|src/**/*.hpp|src/**/*.hxx) return 0 ;;
    *) return 1 ;;
  esac
}

is_excluded_file() {
  case "$1" in
    docs/html/*|third_party/*|build/*|cmake-build-debug/*|cmake-build-release/*|bin/*|obj/*) return 0 ;;
    *.ipynb|*.pdf|*.png|*.jpg|*.jpeg|*.gif|*.svg|*.eps|*.docx) return 0 ;;
    *) return 1 ;;
  esac
}

files=()
if ((all_files)); then
  while IFS= read -r -d '' path; do
    if is_source_file "$path" && ! is_excluded_file "$path"; then
      files+=("$path")
    fi
  done < <(git ls-files -z)
else
  while IFS= read -r -d '' path; do
    if [[ -f "$path" ]] && is_source_file "$path" && ! is_excluded_file "$path"; then
      files+=("$path")
    fi
  done < <(git diff --name-only --diff-filter=ACMR -z "$base_ref"...HEAD)

  while IFS= read -r -d '' path; do
    if [[ -f "$path" ]] && is_source_file "$path" && ! is_excluded_file "$path"; then
      files+=("$path")
    fi
  done < <(git diff --name-only --diff-filter=ACMR -z)

  while IFS= read -r -d '' path; do
    if [[ -f "$path" ]] && is_source_file "$path" && ! is_excluded_file "$path"; then
      files+=("$path")
    fi
  done < <(git ls-files --others --exclude-standard -z)
fi

unique_files=()
if ((${#files[@]} > 0)); then
  for path in "${files[@]}"; do
    seen=0
    for existing in "${unique_files[@]}"; do
      if [[ "$existing" == "$path" ]]; then
        seen=1
        break
      fi
    done
    if ((seen == 0)); then
      unique_files+=("$path")
    fi
  done
fi

if ((${#unique_files[@]} == 0)); then
  echo "No C/C++ files to format."
  exit 0
fi

if [[ "$mode" == "check" ]]; then
  "$clang_format_bin" --dry-run --Werror --style=file "${unique_files[@]}"
else
  "$clang_format_bin" -i --style=file "${unique_files[@]}"
fi

printf '%s\n' "${unique_files[@]}"
