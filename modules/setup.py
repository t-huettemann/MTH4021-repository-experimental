from setuptools import setup

setup(
    name='MTH4021'
    version='0.1.0',
    description="Various Python modules for use with MTH4021",
    author="Thomas Huettemann",
    author_email='t.huettemann@qub.ac.uk',
    python_requires='>=3.7',
    packages = ["rings_and_fields", "abelian_groups", "elliptic_curves"],
    install_requires=['gmpy2', 'primefac']
)
