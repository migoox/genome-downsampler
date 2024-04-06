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

# Install boost library
$USE_SUDO sudo apt install libboost-all-dev