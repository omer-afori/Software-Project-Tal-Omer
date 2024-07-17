from setuptools import Extension, setup

module = Extension("my_capi", sources=['utils.c', 'symnmf.c', 'symnmfmodule.c'])
setup (name="my_capi", version='1.0', description='Python wrapper', ext_modules=[module])