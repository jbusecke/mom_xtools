from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name="mom_xtools",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Collection of xarray based analysis tools for MOM ocean model output",
    author="JuliusBusecke",
    author_email="jbusecke@princeton.edu",
    url="https://github.com/jbusecke/mom_xtools",
    packages=["mom_xtools"],
    entry_points={"console_scripts": ["mom_xtools=mom_xtools.cli:cli"]},
    install_requires=requirements,
    keywords="mom_xtools",
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6",
    ],
)
