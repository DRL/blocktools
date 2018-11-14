from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
ext_modules = [
    Extension("code",  ["code.py"])
    #Extension("mymodule2",  ["mymodule2.py"]),
#   ... all your modules that need be compiled ...
]
setup(
    name = 'blocktools',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)