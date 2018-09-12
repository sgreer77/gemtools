from setuptools import setup 

setup(
        name = 'gemtools',
        packages = ['gemtools'],
        version = '0.1',
        description = 'Tools to analyze 10X data',
        author = 'Hanlee Ji lab',
        author_email = 'sgreer2@stanford.edu',
        url = 'https://github.com/sgreer77/gemtools2',
        dowload_url = 'https://github.com/sgreer77/gemtools/archive/0.1.tar.gz',
        license = 'MIT',
        entry_points={'console_scripts': ['gemtools = gemtools.__main__:main']},
        #package_data = {'SVassembly': ['*.R']},#'*.Rscript']}, #added bc of R file
        #include_package_data=True  #added bc of R file
        #install_requires = ['pandas']  #need to fill this out
        #keywords
        #classifiers
)
