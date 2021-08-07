from distutils.core import setup

setup(
    name='binderspace',
    packages=['binderspace'],
    version='0.1',
    license='mit',
    description='a library for finding the space for binding sequences',
    author='Payam Kelich',
    author_email='payam.kelich68@gmail.com',  # Type in your E-Mail
    url='https://github.com/vukoviclab/BinderSpace',  # Provide either the link to your github or to your website
    download_url='https://github.com/user/reponame/archive/v_01.tar.gz',  # I explain this later on
    keywords=['DNA', 'sequence', 'pca', 't-sne', 'peptide'],  # Keywords that define your package best
    install_requires=[
        'scikit-learn',
        'numpy',
        'pandas',
        'matplotlib',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
)
