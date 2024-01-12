#!/usr/bin/env python

import platform
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler
from setuptools import setup, Extension
from subprocess import getoutput


def using_clang():
    """Will we be using a clang compiler?"""
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler_ver = getoutput("{0} -v".format(compiler.compiler[0]))
    return "clang" in compiler_ver


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked."""

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11

        return pybind11.get_include(self.user)


cpp_extra_link_args = []
cpp_extra_compile_args = [
    "-std=c++11", # C++11
    "-O3",        # O3 optimizations
]

# Needed to specify C++ runtime library on OSX. This solution is replicated
# from the setup.py of mdanalysis
if platform.system() == "Darwin" and using_clang():
    cpp_extra_compile_args.append("-stdlib=libc++")
    cpp_extra_compile_args.append("-mmacosx-version-min=10.9")
    cpp_extra_link_args.append("-stdlib=libc++")
    cpp_extra_link_args.append("-mmacosx-version-min=10.7")

extensions = [
    Extension(
        "matid.ext",
        [
            "matid/ext/ext.cpp",
            "matid/ext/geometry.cpp",
            "matid/ext/celllist.cpp",
        ],
        include_dirs=[
            # Path to pybind11 headers
            "matid/ext",
            get_pybind_include(),
            get_pybind_include(user=True),
        ],
        language="c++",
        extra_compile_args=cpp_extra_compile_args
        + ["-fvisibility=hidden"],  # the -fvisibility flag is needed by pybind11
        extra_link_args=cpp_extra_link_args,
    )
]

if __name__ == "__main__":
    setup(ext_modules=extensions)
