USAGE = "$0 <path to executable> <outputs directory>"

for i in {1..3}; do $1 test -o $2 -t efficiency; done
