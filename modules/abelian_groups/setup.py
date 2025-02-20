from setuptools import setup

setup(
    name='abelian_groups',
    version='0.1.0',
    description="Various Python modules for use with MTH4021",
    author="Thomas Huettemann",
    author_email='t.huettemann@qub.ac.uk',
    python_requires='>=3.7',
    py_modules = ["abelian_groups"],
    install_requires=['gmpy2']
)
