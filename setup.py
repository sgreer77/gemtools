from setuptools import setup 

setup(
        name = 'gemtools2',
        packages = ['gemtools2'],
        version = '0.2',
        description = 'Tools to analyze 10X data',
        author = 'Hanlee Ji lab',
        author_email = 'sgreer2@stanford.edu',
        url = 'https://github.com/sgreer77/gemtools2',
        dowload_url = 'https://github.com/sgreer77/gemtools2/archive/0.2.tar.gz',
        license = 'MIT',
        entry_points={'console_scripts': ['gemtools2 = gemtools2.__main__:main']},
        #package_data = {'SVassembly': ['*.R']},#'*.Rscript']}, #added bc of R file
        #include_package_data=True  #added bc of R file
        #install_requires = ['pandas']  #need to fill this out
        #keywords
        #classifiers
)
