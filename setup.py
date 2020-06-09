import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="starpolymers",
    version="2.0.0",
    author="Debesh Mandal",
    description="Package for creating polymers for simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/debeshmandal/starpolymers",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.5",
)
