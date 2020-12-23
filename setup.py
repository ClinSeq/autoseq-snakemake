from setuptools import setup, find_packages

setup(name='autoseq',
      version='0.1.0',
      packages=find_packages(exclude=('tests*', 'docs', 'examples')),
      entry_points={
          'console_scripts': [
              'autoseq = pipeline.cli:cli',
          ]
      }
    )
