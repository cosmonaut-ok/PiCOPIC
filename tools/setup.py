from setuptools import setup
# from distutils.core import setup

setup(name='pdp3-tools',
      version='0.1.0',
      description='Data analysis tools for PDP3 project',
      author='Alexander Vynnyk',
      author_email='alexander.vynnyk@gmail.com',
      url='https://gitlab.com/my-funny-plasma/PIC/pdp3',
      install_requires=[
          'matplotlib', 'numpy', 'scipy', 'h5py'
      ],
      zip_safe=False)
