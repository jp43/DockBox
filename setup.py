import sys
import setuptools

def exit_with_error(head, body=''):
    _print_admonition('error', head, body)
    sys.exit(1)

# check Python version
if not (sys.version_info[0] == 2 and sys.version_info[1] >= 6):
    exit_with_error("You need Python 2.6.x or Python 2.7.x to install the DockBox package!")

setuptools.setup(name='dockbox',
    version='1.1',
    packages=['dockbox'],
    scripts=['bin/rundbx', 'bin/extract_dbx_best_poses'],
    install_requires=['mdkit', 'cython', 'numpy==1.8.0', 'pandas==0.19.0'],
    license='LICENSE.txt',
    description='Platform package to simplify the use of docking programs and consensus methods',
    long_description=open('README.rst').read(),
)
