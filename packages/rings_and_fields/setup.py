from setuptools import setup

setup(
    name='rings_and_fields',
    version='0.1.0',
    description="Rings and Fields module",
    author="Thomas Huettemann",
    author_email='t.huettemann@qub.ac.uk',
    url="https://github.com/t-huettemann/MTH4021-repository/packages/rings_and_fields",
    python_requires='>=3.7',
    modules=["rings_and_fields"]
    #packages=["rings_and_fields", "abelian_groups", "elliptic_curves"],
    #package_dir={"rings_and_fields":"packages/rings_and_fields", "abelian_groups":"packages/abelian_groups", "elliptic_curves":"packages/elliptic_curves"},
    install_requires=['gmpy2', 'primefac']
)
