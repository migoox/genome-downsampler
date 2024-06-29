USAGE = "$0 <path to executable> <outputs directory>"

# for i in {1..10}; do $1 test -o $2 -t efficiency -a quasi-mcp-cpu; done
for i in {1..10}; do time $1 test -o $2 -t efficiency -a qmcp-cpu mcp-cpu quasi-mcp-cuda; done

