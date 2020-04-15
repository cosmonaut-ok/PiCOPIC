from setuptools import setup
# from distutils.core import setup

setup(name='picopic',
      version='0.1.0',
      description='Data analysis tools for PiCOPIC project',
      author='Alexander Vynnyk',
      author_email='alexander.vynnyk@gmail.com',
      url='https://gitlab.com/my-funny-plasma/PIC/picopic',
      install_requires=[
          'matplotlib', 'numpy', 'scipy', 'h5py', 'setupext-janitor', 'pygments'
      ],
      zip_safe=False)
