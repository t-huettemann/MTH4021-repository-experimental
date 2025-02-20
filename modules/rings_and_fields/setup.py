from setuptools import setup

setup(
    name='rings_and_fields',
    version='0.1.0',
    description="Various Python modules for use with MTH4021",
    author="Thomas Huettemann",
    author_email='t.huettemann@qub.ac.uk',
    python_requires='>=3.7',
    py_modules = ["rings_and_fields"],
    install_requires=['gmpy2', 'primefac']
)
