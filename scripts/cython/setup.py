from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension("edgelistParser",  ["edgelistParser.pyx"]),
    Extension("pageRank",  ["pageRank.pyx"]),
    Extension("plotNetwork",  ["plotNetwork.pyx"]),
    Extension("utils",  ["utils.pyx"])
]

setup(
    name = 'LocalPagerank Utils',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
