from setuptools import setup, find_packages

setup(
    name="blastoise",
    version="0.4.2",
    author="R. Pacheco",
    description="A Tool for Repetitive Sequence Discovery in Genomic Data",
    packages=find_packages(),
    python_requires=">=3.10",
    # No install_requires - let conda handle all dependencies
    entry_points={
        "console_scripts": [
            "blastoise=blastoise.main:main",
            "blastoise_coordinate_corrector=blastoise.extra.coordinate_corrector:main",
            "blastoise_sider_filter=blastoise.extra.sider_filter:main",

        ],
    },
)