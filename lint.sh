# git diff -U0 HEAD^ | ./scripts/diff-clang-tidy.py -p1 -path ./build/
git diff -U0 HEAD^ | ./scripts/run-clang-tidy.py -p=./build/