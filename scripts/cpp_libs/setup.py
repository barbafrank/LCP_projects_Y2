from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

#from Cython.Build import cythonize

ext_modules = [
    Extension("pageRankCpp",  ["pageRankCpp.pyx"], extra_compile_args = ["-stdlib=libc++","-std=c++11","-Ofast","-march=native", "-mmacosx-version-min=10.9"], language="c++", extra_link_args=["-stdlib=libc++", "-mmacosx-version-min=10.9"])
]

setup(
    name = 'LocalPagerank Utils',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
