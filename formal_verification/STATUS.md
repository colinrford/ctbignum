# Formal Verification Status (December 2025)
# Disclaimer: This document wasargely written by Gemini

This directory contains the formal verification wrappers and scripts for **ctbignum**.

## Summary
The verification code in this repository has been fully modernized to support **C++20 Modules** and **C++23** standards. However, as of Dec 2025, the external verification tools (**SAW** and **ct-verif**) are not yet compatible with the **LLVM 19** bitcode required by this modern build.

## Code Status: READY
All C++ wrappers have been updated to use C++20 module syntax (`import ctbignum;`, `import std;`) and compile successfully to LLVM bitcode using `clang++-19`.

-   **SAW Wrappers** (`saw-cryptol/`):
    -   `add.cpp`, `mul.cpp`, `mul_wrapper.cpp`: **Updated & Compiles**
-   **ct-verif Wrappers** (`ct-verif/`):
    -   `wrapper.cpp`: **Updated & Compiles**

## Tool Status: BLOCKED
The current latest versions of the verification tools (e.g., SAW v1.4) crash when parsing LLVM 19 bitcode.

-   **Error**: `FUNC_CODE_INST_GEP` (Bitcode parsing error)
-   **Cause**: Clang 19 generates a newer LLVM bitcode format that the Haskell-based LLVM parsers in these tools do not yet support.

## How to Run (Future)
Once the tools are updated to support LLVM 19:

1.  **Update Dockerfile**: Change the SAW/ct-verif installation to the new compatible version.
2.  **Run Script**:
    ```bash
    docker run --rm -v $(pwd):/app -w /app ctbignum-modern bash formal_verification/saw-cryptol/verify.sh
    ```
