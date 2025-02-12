from setuptools import setup, find_packages

setup(
    name='CFG',
    version='1.0.0',
    author='',
    author_email='',
    description='AI calculater using smallest grammar algorithm re-pair.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/ELIFE-ASU/CFG',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12',
    install_requires=[
        'networkx>=3.4.2'
    ],
    extras_require={
        'dev': [
            'pytest',
            'pytest-cov',
        ],
    },
)
