from skbuild import setup
from setuptools import find_packages
from pathlib import Path
import os

# print(os.environ["CMAKE_PREFIX_PATH"])
p = "Plugins/Python/python"
# p = "python"
packages = find_packages(where=p)
print(packages)
print(p)
# import sys
# sys.exit()

cmake_args = [
    "-DACTS_BUILD_EXAMPLES=ON",
    "-DCMAKE_CXX_STANDARD=17",
    "-DCMAKE_CXX_FLAGS=-fdiagnostics-color=always -fPIC",
    "-DACTS_BUILD_PLUGIN_PYTHON=ON",
    "-DCMAKE_BUILD_TYPE=Debug",
]

if "CMAKE_PREFIX_PATH" in os.environ:
    cmake_args.append("-DCMAKE_PREFIX_PATH=" + os.environ["CMAKE_PREFIX_PATH"])

setup(
    name="acts",
    version="0.1",
    description="ACTS",
    cmake_args=cmake_args,
    packages=["acts"],
    package_dir={"": p},
    cmake_install_dir="Plugins/Python/python/acts",
    #     install_requires=[],
    #     package_data = {
    #         'acts': ["*.so"]
    #     },
    include_package_data=True,
)
