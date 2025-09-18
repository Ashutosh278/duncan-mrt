from setuptools import setup, find_packages

# Read the README.md file correctly with UTF-8 encoding
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="duncans-mrt",
    version="0.1.0",
    author="Ashutosh Sahu",
    description="Duncan's Multiple Range Test implementation with plotting",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
    ],
    keywords="statistics, multiple comparison, duncan test, anova",
)
