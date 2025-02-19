from setuptools import setup

setup(
    name='rings_and_fields',
    version='1.0',
    description="Rings and Fields package",
    author="Thomas Huettemann",
    author_email='t.huettemann@qub.ac.uk',
    url="https://github.com/t-huettemann/MTH4021-repository/",
    python_requires='>=3.7',
    packages=["rings_and_fields"],
    package_dir={"rings_and_fields":"packages/rings_and_fields"},
    install_requires=['gmpy2', 'primefac']
)
