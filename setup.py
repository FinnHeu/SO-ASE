from setuptools import setup, find_packages

setup(
    name="so-ase",
    version="0.1.0",
    description="Standard analysis toolbox.",
    author="Finn Ole Heukamp",
    author_email="finn.heukamp@awi.de",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'cartopy',
        'xarray',
        'pyfesom2',
        'scipy',
        'ipykernel'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.7",
)
