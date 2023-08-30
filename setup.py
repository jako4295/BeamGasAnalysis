"""
setup.py for project-start.

For reference see
https://packaging.python.org/guides/distributing-packages-using-setuptools/

"""
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
VERSION = "0.0.1"
LONG_DESCRIPTION = (here / "README.md").read_text(encoding="utf-8")

REQUIREMENTS: dict = {
    "core": [
        "numpy",
        "matplotlib",
        "scipy",
        "pandas",
    ],
    "test": [],
    "dev": [
        # 'requirement-for-development-purposes-only',
    ],
    "doc": [],
}


setup(
    name="BeamGas",
    version=VERSION,
    author="Jakob Olsen",
    author_email="jakob.olsen@cern.ch",
    description="Calculation of beam gas lifetimes and cross sections.",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="",
    packages=find_packages(),
    python_requires="~=3.11.0",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=REQUIREMENTS["core"],
    extras_require={
        **REQUIREMENTS,
        # The 'dev' extra is the union of 'test' and 'doc', with an option
        # to have explicit development dependencies listed.
        "dev": [
            req
            for extra in ["dev", "test", "doc"]
            for req in REQUIREMENTS.get(extra, [])
        ],
        # The 'all' extra is the union of all requirements.
        "all": [req for reqs in REQUIREMENTS.values() for req in reqs],
    },
    include_package_data=True,
)
