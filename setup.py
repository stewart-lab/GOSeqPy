from setuptools import setup

setup(
    name = 'GOSeqPy', 
    version = '0.0.1',
    description = (
        'A python wrapper for running the GOSeq R Package'
    ),
    py_modules = ["run_goseq"],
    package_dir = {'':'src'},
    author = 'Ron Stewart',
    author_email = 'rstewart@morgridge.org',
    long_description = open('README.md').read() + '\n\n' + open('CHANGELOG.md').read(),
    long_description_content_type = "text/markdown",
    url='https://github.com/stewart-lab/GOSeqPy',
    include_package_data=True,
    install_requires = ['numpy', 'cmdlogtime', 'pandas', 'cmdlogtime'],
    keywords = ['GOSeq'],
)

# Whenever you want to update your package, you should remove the ‘build’ and
# ‘dist’ folders, make changes to your code, edit the “CHANGLOG.txt” file, and 
# revise the version number in the “setup.py”.
#
# conda activate cmdlogtime
# rm -rf build/ dist/
# python setup.py sdist bdist_wheel
# pytest (NA yet)
# twine check dist/*
# twine upload --repository-url https://test.pypi.org/legacy/ dist/* 
#   #pay attention there is an extra space before dist.
# twine upload dist/*
