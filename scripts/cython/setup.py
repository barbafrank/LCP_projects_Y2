from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [
    Extension("edgelistParser",  ["edgelistParser.pyx"], extra_compile_args = ["-Ofast -march=native"]),
    Extension("pageRank",  ["pageRank.pyx"], include_dirs=[numpy.get_include()], extra_compile_args = ["-Ofast -march=native"]),
    Extension("plotNetwork",  ["plotNetwork.pyx"], extra_compile_args = ["-Ofast -march=native"]),
    Extension("utils",  ["utils.pyx"], include_dirs=[numpy.get_include()], extra_compile_args = ["-Ofast -march=native"])
]

setup(
    name = 'LocalPagerank Utils',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
