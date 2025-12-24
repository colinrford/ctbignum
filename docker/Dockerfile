FROM ubuntu:24.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install basics and ninja
RUN apt-get update && apt-get install -y \
    wget \
    git \
    ninja-build \
    software-properties-common \
    lsb-release \
    gnupg \
    libc++-19-dev \
    libc++abi-19-dev \
    z3 \
    libgmp-dev \
    libntl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install LLVM 19 (Support for C++23 modules)
RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)" -- 19
ENV CC=clang-19
ENV CXX=clang++-19

# Install CMake 3.31.6 (Required by CMakeLists.txt)
RUN wget https://github.com/Kitware/CMake/releases/download/v3.31.6/cmake-3.31.6-linux-x86_64.sh \
    && chmod +x cmake-3.31.6-linux-x86_64.sh \
    && ./cmake-3.31.6-linux-x86_64.sh --prefix=/usr/local --exclude-subdir --skip-license \
    && rm cmake-3.31.6-linux-x86_64.sh

# Install SAW (Software Analysis Workbench) v1.4
RUN wget https://github.com/GaloisInc/saw-script/releases/download/v1.4/saw-1.4-ubuntu-24.04-X64-with-solvers.tar.gz \
    && tar -xf saw-1.4-ubuntu-24.04-X64-with-solvers.tar.gz \
    && mv saw-1.4-ubuntu-24.04-X64-with-solvers/bin/* /usr/local/bin/ \
    && mv saw-1.4-ubuntu-24.04-X64-with-solvers/lib/* /usr/local/lib/ \
    && rm -rf saw-1.4-ubuntu-24.04-X64-with-solvers*

# Build Google Benchmark from source (libc++ compatible)
RUN git clone https://github.com/google/benchmark.git \
    && cd benchmark \
    && cmake -E make_directory "build" \
    && cmake -E chdir "build" cmake -GNinja -DBENCHMARK_DOWNLOAD_DEPENDENCIES=on -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++-19 -DCMAKE_CXX_FLAGS="-stdlib=libc++" -DCMAKE_EXE_LINKER_FLAGS="-stdlib=libc++ -lc++abi" ../ \
    && cmake --build "build" --config Release \
    && cmake --install "build"

WORKDIR /app
COPY . .

# Verification command (can be overridden)
CMD ["bash", "-c", "rm -rf build && cmake -G Ninja -B build -DCTBIGNUM_BuildTests=ON -DCMAKE_CXX_FLAGS=\"-stdlib=libc++\" -DCMAKE_EXE_LINKER_FLAGS=\"-stdlib=libc++\" . && cmake --build build && cd build && ctest"]
