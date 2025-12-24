#!/bin/bash
set -e

echo "=== 1. Building Library (Modules) ==="
rm -rf build
cmake -G Ninja -B build -DCTBIGNUM_BuildTests=ON -DCMAKE_CXX_FLAGS="-stdlib=libc++" -DCMAKE_EXE_LINKER_FLAGS="-stdlib=libc++" .
cmake --build build

echo "=== 2. Compiling Verification Wrappers to Bitcode ==="
FLAGS="-std=c++23 -stdlib=libc++ -fprebuilt-module-path=build/CMakeFiles/ctbignum.dir -fprebuilt-module-path=build/CMakeFiles/__cmake_cxx23.dir -c -emit-llvm"

echo "Compiling add.cpp..."
clang++-19 $FLAGS formal_verification/saw-cryptol/add.cpp -o formal_verification/saw-cryptol/add.bc

echo "Compiling mul.cpp..."
clang++-19 $FLAGS formal_verification/saw-cryptol/mul.cpp -o formal_verification/saw-cryptol/mul.bc

echo "Compiling mul_wrapper.cpp..."
clang++-19 $FLAGS formal_verification/saw-cryptol/mul_wrapper.cpp -o formal_verification/saw-cryptol/mul_wrapper.bc

echo "=== 3. Running Formal Verification (SAW) ==="
cd formal_verification/saw-cryptol

echo "Verifying Addition..."
saw add_crucible.saw

echo "Verifying Muliplication (Crucible)..."
saw mul_crucible.saw

echo "Verifying Multiplication (Wrapper)..."
saw mul.saw

echo "=== Verification Complete! ==="
