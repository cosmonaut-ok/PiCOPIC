from setuptools import setup
# from distutils.core import setup

setup(name='picopic-test',
      version='0.1.0',
      description='Functional Testing framework for PiCoPiC project',
      author='Alexander Vynnyk',
      author_email='alexander.vynnyk@gmail.com',
      url='https://gitlab.com/my-funny-plasma/PIC/picopic',
      install_requires=[
          'numpy', 'h5py', 'pygments', 'jinja2', 'colorama'
      ],
      zip_safe=False)
