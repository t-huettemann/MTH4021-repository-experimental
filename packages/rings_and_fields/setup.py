from setuptools import setup, find_packages

setup(
    name='rings_and_fields',
    version='1.0',
    description="Rings and Fields package",
    author="Thomas Huettemann",
    author_email='t.huettemann@qub.ac.uk',
    url="https://github.com/t-huettemann/MTH4021-repository.git",
    python_requires='>=3.7',
    packages=find_packages(),
    install_requires=['gmpy2', 'primefac']
)
