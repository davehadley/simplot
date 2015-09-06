from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

language = "c++"
extra_compile_args=["-std=c++11", "-O3"]
extra_link_args=["-std=c++11", "-O3"]


modules = [Extension("simplot.rootplot.cmerge",
                     ["simplot/rootplot/cmerge.pyx"],
                     language = language,
                     extra_compile_args=extra_compile_args,
                     extra_link_args=extra_link_args),
#           Extension("simplot.mc.generators",
#                     ["simplot/mc/generators.pyx"],
#                     language = language,
#                     extra_compile_args=extra_compile_args,
#                     extra_link_args=extra_link_args),
           Extension("simplot.mc.clikelihood",
                     ["simplot/mc/clikelihood.pyx"],
                     language = language,
                     extra_compile_args=extra_compile_args,
                     extra_link_args=extra_link_args),
           Extension("simplot.sparsehist.sparsehist",
                     ["simplot/sparsehist/unordered_map.pxd", "simplot/sparsehist/sparsehist.pxd", "simplot/sparsehist/sparsehist.pyx"],
                     language = language,
                     extra_compile_args=extra_compile_args,
                     extra_link_args=extra_link_args),
           Extension("simplot.binnedmodel.model",
                     ["simplot/binnedmodel/model.pyx"],
                     language = language,
                     extra_compile_args=extra_compile_args,
                     extra_link_args=extra_link_args),
           Extension("simplot.binnedmodel.fluxweights",
                 ["simplot/binnedmodel/fluxweights.pyx"],
                 language = "c++",
                     extra_compile_args=["-std=c++11", "-O3"],
                     extra_link_args=["-std=c++11", "-O3"]),
           Extension("simplot.binnedmodel.xsecweights",
                 ["simplot/binnedmodel/xsecweights.pyx"],
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
