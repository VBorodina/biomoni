import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name= "biomoni",
    version='0.0.1',
    author= "Paul",
    packages=["biomoni"],
    license="LICENSE.txt",
    description="Package to automate Biomonitoring",
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        "Main_repo": "https://github.com/PSenck/Biomonitoring",
    }
)