from setuptools import setup 

setup(
        name = 'gemtools',
        packages = ['gemtools'],
        version = '1.0',
        description = 'Tools to analyze 10X data',
        author = 'Stephanie Greer / Ji Research Group',
        author_email = 'sgreer2@stanford.edu',
        url = 'https://github.com/sgreer77/gemtools',
        dowload_url = 'https://github.com/sgreer77/gemtools/archive/1.0.tar.gz',
        license = 'MIT',
        entry_points={'console_scripts': ['gemtools = gemtools.__main__:main']},
)
