#!/bin/zsh

echo "=== Building config files for debug mode ==="
cmake -S . -B build/debug -G Ninja -DCMAKE_BUILD_TYPE=debug

echo "\n=== Building config files for release mode ==="
cmake -S . -B build/release -G Ninja -DCMAKE_BUILD_TYPE=release

