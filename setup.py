from setuptools import setup, find_packages

version = '0.2.0'

setup(name='autoseq',
      version=version,
      packages=find_packages(exclude=('tests*', 'docs', 'examples')),
      install_requires=[
          "snakemake",
          "click",
          "pyyaml",
          "pandas",
          "rich",
          "loguru"
      ],
      entry_points={
          'console_scripts': [
              'autoseq = pipeline.cli:console_autoseq',
          ]
      }
    )
