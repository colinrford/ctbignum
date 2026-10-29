# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2018 Niek J. Bouman (original CTBignum)
# SPDX-FileCopyrightText: 2025-2026 Colin Ford (lam modularization)
#
# Conan 2.x recipe for lam.ctbignum: conan create . --profile <profile>

from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout
from conan.tools.files import copy, load
import os

class LamCtbignumConan(ConanFile):
    name = "lam_ctbignum"
    license = "Apache-2.0"
    author = "Colin Ford"
    url = "https://github.com/colinrford/ctbignum"
    description = (
        "Compile-time and run-time multi-precision and modular arithmetic in "
        "C++23 modules for the lam project. A modularized fork of Niek Bouman's "
        "CTBignum."
    )
    topics = ("c++23", "modules", "bignum", "modular-arithmetic", "constexpr",
              "cryptography", "lam")

    settings = "os", "compiler", "build_type", "arch"
    package_type = "static-library"

    def set_version(self):
        self.version = load(
            self, os.path.join(self.recipe_folder, "VERSION")
        ).strip()

    exports_sources = (
        "VERSION",
        "CMakeLists.txt",
        "include/*",
        "cmake/*",
        "LICENSE",
        "NOTICE",
        "README.md",
    )

    def layout(self):
        cmake_layout(self)

    def validate(self):
        cppstd = self.settings.compiler.cppstd
        if cppstd is not None:
            std = int(str(cppstd).replace("gnu", ""))
            if std < 23:
                raise Exception(
                    "lam_ctbignum requires C++23 (compiler.cppstd >= 23)."
                )

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["CMAKE_CXX_STANDARD"] = "23"
        tc.cache_variables["CMAKE_CXX_SCAN_FOR_MODULES"] = "ON"
        # test/ and benchmarks/ aren't exported
        tc.cache_variables["LAM_CTBIGNUM_BUILD_TESTS"] = "OFF"
        tc.cache_variables["LAM_CTBIGNUM_BUILD_BENCHMARKS"] = "OFF"
        tc.generate()
        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()
        for f in ("LICENSE", "NOTICE"):
            copy(
                self,
                f,
                src=self.source_folder,
                dst=os.path.join(self.package_folder, "licenses"),
            )

    def package_info(self):
        self.cpp_info.set_property("cmake_file_name", "lam_ctbignum")
        self.cpp_info.set_property("cmake_target_name", "lam_ctbignum::ctbignum")
        self.cpp_info.libs = ["lam_ctbignum"]
        self.cpp_info.builddirs = [
            os.path.join("lib", "cmake", "lam_ctbignum"),
            os.path.join("lib", "cps"),
        ]
