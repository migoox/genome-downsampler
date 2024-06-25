# Detect if the script is being run as root, storing true/false in is_root.
is_root=false
if (( $EUID == 0)); then
   is_root=true
fi
# Find if sudo is available
has_sudo=false
if command -v sudo &> /dev/null ; then
    has_sudo=true
fi
# Decide if we can proceed or not (root or sudo is required) and if so store whether sudo should be used or not. 
if [ "$is_root" = false ] && [ "$has_sudo" = false ]; then 
    echo "Root or sudo is required. Aborting."
    exit 1
elif [ "$is_root" = false ] ; then
    USE_SUDO=sudo
else
    USE_SUDO=
fi


$USE_SUDO apt update 

# Install htslib dependencies (https://github.com/samtools/htslib/blob/develop/INSTALL)
$USE_SUDO apt install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

# Install htslib only
# By default, 'make install' installs HTSlib libraries under /usr/local/lib,
# HTSlib header files under /usr/local/include
(
    git clone --recurse-submodules https://github.com/samtools/htslib.git
    cd htslib
    make
    $USE_SUDO make install
)

# Install samtools and htslib using the package manager
# $USE_SUDO apt install samtools

# Install ortools under /usr/local/
ORTOOLS_VERSION=or-tools_amd64_ubuntu-22.04_cpp_v9.9.3963
ORTOOLS_DIR_NAME=or-tools_x86_64_Ubuntu-22.04_cpp_v9.9.3963
(
    wget https://github.com/google/or-tools/releases/download/v9.9/${ORTOOLS_VERSION}.tar.gz
    tar -xvzf ${ORTOOLS_VERSION}.tar.gz
    $USE_SUDO cp -r ${ORTOOLS_DIR_NAME}/bin/* /usr/local/bin/
    $USE_SUDO cp -r ${ORTOOLS_DIR_NAME}/lib/* /usr/local/lib/
    $USE_SUDO cp -r ${ORTOOLS_DIR_NAME}/include/* /usr/local/include/
    $USE_SUDO cp -r ${ORTOOLS_DIR_NAME}/share/* /usr/local/share/
)