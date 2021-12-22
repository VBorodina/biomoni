# How to install the biomoni library?

Create new anaconda environment (optional for other env than conda base):

`conda create –name Biomonitoring python`

Activate environment

`conda activate Biomonitoring`

Copy the folder to your local dirve

To install the required packages: open directory in terminal and type:

`pip install -r Requirements.txt`

To install the package in the current Anaconda env:git add **

`pip install -e .`      # . means this directory

To install from github directly without cloning the folder:

Requirements:
`pip install -r https://raw.githubusercontent.com/PSenck/biomoni/main/Requirements.txt`

biomoni:
`pip install -e git+https://github.com/PSenck/biomoni.git#egg=biomoni`





