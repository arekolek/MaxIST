
from distutils.core import setup
from distutils.extension import Extension

setup(
  name="maxist",
  version="0.1",
  author="Arkadiusz Olek",
  author_email="arkadiusz.olek@uj.edu.pl",
  ext_modules=[
    Extension(
      name = "graph",
      sources = ["graph.cpp"],
      include_dirs = ["graph", "util", "test"],
      extra_compile_args = ["-std=c++0x"],
      libraries = ["boost_python-py34"])
  ])
