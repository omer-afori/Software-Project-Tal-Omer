from setuptools import Extension, setup

module = Extension("symnmfmodule", sources=['utils.c', 'symnmf.c', 'symnmfmodule.c'])
setup (name="symnmfmodule", version='1.0', description='Python wrapper', ext_modules=[module])