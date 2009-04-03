from distutils.core import setup, Extension

setup(name='SnppnS',
      version='\omega+1',
      ext_modules=[Extension('_Dirichlet', ['_Dirichlet.c'])],
      )
