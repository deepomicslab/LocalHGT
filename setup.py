from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext


from Cython.Compiler.Options import get_directive_defaults
directive_defaults = get_directive_defaults()

directive_defaults['linetrace'] = True
directive_defaults['binding'] = True

#setup(
#    ext_modules=cythonize("speed_test.pyx"),
#)


ext_modules=[ Extension("sample_specific_genomes_v2",
              ["sample_specific_genomes_v2.pyx"],
              libraries=["m"],
              extra_compile_args=["-ffast-math"])]
              
setup(
    name='sample_specific_genomes_v2', 
    cmdclass={"build_ext": build_ext},
    ext_modules=ext_modules)


extensions = [
    Extension("sample_specific_genomes_v2", ["sample_specific_genomes_v2.pyx"], define_macros=[('CYTHON_TRACE', '1')])
]
