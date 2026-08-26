from setuptools import setup, find_packages

version = '4.1.3'

setup(name='autoseq',
      version=version,
      packages=find_packages(exclude=('tests*', 'docs', 'examples')),
      python_requires=">=3.11",
      install_requires=[
          "snakemake==9.23.0",
          "snakemake-executor-plugin-slurm==2.7.1",
          "click",
          "pyyaml",
          "pandas",
          "rich",
          "loguru",
          "pytest",
      ],
      entry_points={
          'console_scripts': [
              'autoseq = pipeline.cli:console_autoseq',
          ]
      }
    )
