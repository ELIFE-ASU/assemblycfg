from setuptools import setup, find_packages

setup(
    name='CFG',
    version='0.1.0',
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
    python_requires='>=3.9',
    install_requires=[
    ],
    extras_require={
        'dev': [
            'pytest',
            'pytest-cov',
        ],
    },
)
