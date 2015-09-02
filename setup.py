from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"


modules = [Extension("simplot.rootplot.cmerge",
                     ["simplot/rootplot/cmerge.pyx"],
                     language = "c++",
                     extra_compile_args=["-std=c++11", "-O3"],
                     extra_link_args=["-std=c++11", "-O3"]),
]

setup(name="simplot",
      version="0.1",
      description="Set of tools for plotting with matplotlib and ROOT.",
      author="David Hadley",
      author_email="d.r.hadley@warwick.ac.uk",
      url="https://github.com/davehadley/simplot",
      packages=["simplot", "simplot.batch", "simplot.mplot", "simplot.rootplot"],
      cmdclass={"build_ext": build_ext},
      ext_modules=modules,
)
