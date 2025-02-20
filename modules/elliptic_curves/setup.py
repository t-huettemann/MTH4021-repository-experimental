from setuptools import setup

setup(
    name='elliptic_curves',
    version='0.1.0',
    description="Various Python modules for use with MTH4021",
    author="Thomas Huettemann",
    author_email='t.huettemann@qub.ac.uk',
    python_requires='>=3.7',
    py_modules = ["elliptic_curves"],
    install_requires=['abelian_groups']
)
