from setuptools import setup

setup(
    packages=['adnator'],
    package_dir={'adnator': 'src/adnator'},
    package_data={'adnator': ['*.so']},
    include_package_data=True
)
