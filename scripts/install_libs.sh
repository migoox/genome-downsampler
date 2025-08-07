#!/usr/bin/env bash

set -euo pipefail

DEFAULT_INSTALL_DIR="/usr/local"

function print_help() {
    cat << EOF
Usage: $0 install [--prefix <install_dir>]

Downloads and installs Google OR-Tools and HTSlib C++ binaries based on your Linux distro and CPU architecture.
The script uses the following releases 
    - https://github.com/google/or-tools/releases/tag/v9.9
    - https://github.com/samtools/htslib/releases/tag/1.22.1

Requirements: gcc, make, zlib, libbz2, liblzma, libcurl 

Options:
  install               Download and install OR-Tools binary
  --prefix <dir>        Optional installation directory (default: ${DEFAULT_INSTALL_DIR})
  --subset <components>...    Install only specified components (ortools, htslib). Default: all
  --help                Show this help message

Examples:
  $0 install
  $0 install --prefix /opt/or-tools
EOF
}

function get_cpu_arch() {
    local arch
    arch="$(uname -m | tr '[:upper:]' '[:lower:]')"
    case "$arch" in
        x86_64)
            echo "amd64"
            ;;
        aarch64)
            echo "aarch64"
            ;;
        *)
            echo "Unsupported architecture: $arch" >&2
            exit 1
            ;;
    esac
}

function get_linux_distro() {
    if [[ ! -f /etc/os-release ]]; then
        echo "Cannot detect Linux distro: /etc/os-release not found." >&2
        exit 1
    fi

    source /etc/os-release
    local id="${ID,,}"
    local version_id="${VERSION_ID:-unknown}"

    case "$id" in
        ubuntu)
            echo "ubuntu-${version_id}"
            ;;
        debian)
            if [[ "$version_id" == "unstable" ]]; then
                echo "debian-sid"
            else
                echo "debian-${version_id}"
            fi
            ;;
        almalinux)
            echo "almalinux-${version_id}"
            ;;
        rocky)
            echo "rockylinux-${version_id}"
            ;;
        fedora)
            echo "fedora-${version_id}"
            ;;
        opensuse-leap)
            echo "opensuse-leap"
            ;;
        alpine)
            echo "alpine-edge"
            ;;
        arch | endeavouros | manjaro)
            echo "archlinux"
            ;;
        nixos)
            echo "WARNING: NixOS is not FHS-compliant. OR-Tools binaries likely won't work." >&2
            echo "nixos-unsupported"
            ;;
        *)
            echo "Unsupported Linux distro: $id" >&2
            exit 1
            ;;
    esac
}

function install_ortools() {
    local install_dir="$1"

    echo "[*] Detecting architecture and Linux distribution..."
    local arch
    arch="$(get_cpu_arch)"
    local distro
    distro="$(get_linux_distro)"

    echo "[*] Detected architecture: $arch"
    echo "[*] Detected distro: $distro"

    local version="v9.9"
    local tag="or-tools_${arch}_${distro}_cpp_v9.9.3963"
    local url="https://github.com/google/or-tools/releases/download/${version}/${tag}.tar.gz"

    echo "[*] Downloading OR-Tools from:"
    echo "    $url"

    tmpdir="$(mktemp -d)"
    trap 'rm -rf "$tmpdir"' EXIT

    curl -L "$url" -o "$tmpdir/ortools.tar.gz"
    mkdir -p "$tmpdir/extract"
    ortools_root_dir="$(tar tzf "$tmpdir/ortools.tar.gz" | sed -e 's@/.*@@' | uniq)"
    tar -xf "$tmpdir/ortools.tar.gz" -C "$tmpdir/extract"

    echo "[*] Installing OR-Tools to: $install_dir"
    sudo mkdir -p "$install_dir"
    sudo mkdir -p "$install_dir/lib"
    sudo mkdir -p "$install_dir/include"
    sudo mkdir -p "$install_dir/bin"
    sudo mkdir -p "$install_dir/share"
    sudo cp -r "$tmpdir/extract/$ortools_root_dir/"lib/* "$install_dir/lib"
    sudo cp -r "$tmpdir/extract/$ortools_root_dir/"include/* "$install_dir/include"
    sudo cp -r "$tmpdir/extract/$ortools_root_dir/"bin/* "$install_dir/bin"
    sudo cp -r "$tmpdir/extract/$ortools_root_dir/"share/* "$install_dir/share"

    echo "[✓] OR-Tools installed successfully to $install_dir"
}

function install_htslib() {
    local install_dir="$1"
    local url="https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2"

    echo "[*] Downloading HTSlib from:"
    echo "    $url"

    tmpdir="$(mktemp -d)"
    trap 'rm -rf "$tmpdir"' EXIT

    curl -L "$url" -o "$tmpdir/htslib.tar.bz2"
    mkdir -p "$tmpdir/extract"
    tar -xf "$tmpdir/htslib.tar.bz2" -C "$tmpdir/extract"
    htslib_root_dir="$(tar tjf "$tmpdir/htslib.tar.bz2" | sed -e 's@/.*@@' | uniq)"
    (
        cd "$tmpdir/extract/$htslib_root_dir"
        echo "[*] Compiling HTSlib"
        make

        echo "[*] Installing HTSlib to: $install_dir"
        sudo make DESTDIR=$install_dir prefix="" install 
    )

    echo "[✓] HTSlib installed successfully to $install_dir"
}

# === Entry Point ===

if [[ $# -eq 0 || "$1" == "--help" ]]; then
    print_help
    exit 0
fi

if [[ "$1" == "install" ]]; then
    shift
    install_path="$DEFAULT_INSTALL_DIR"
    subset=()  # empty means install all

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --prefix)
                if [[ -n "${2:-}" ]]; then
                    install_path="$2"
                    shift 2
                else
                    echo "Error: --prefix requires a path argument" >&2
                    exit 1
                fi
                ;;
            --subset)
                shift
                while [[ $# -gt 0 && "$1" != --* ]]; do
                    subset+=("$1")
                    shift
                done
                ;;
            *)
                echo "Unknown argument: $1" >&2
                print_help
                exit 1
                ;;
        esac
    done

    # Install all if subset is empty
    if [[ ${#subset[@]} -eq 0 ]]; then
        subset=("ortools" "htslib")
    fi

    for component in "${subset[@]}"; do
        case "$component" in
            ortools)
                install_ortools "$install_path"
                ;;
            htslib)
                install_htslib "$install_path"
                ;;
            *)
                echo "Unknown subset component: $component" >&2
                exit 1
                ;;
        esac
    done

else
    echo "Unknown command: $1" >&2
    print_help
    exit 1
fi
