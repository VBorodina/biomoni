# How to install the biomoni library?

Create new anaconda environment (optional for other env than conda base):

`conda create –name Biomonitoring python`

Activate environment

`conda activate Biomonitoring`

Copy the folder to your local drive

To install the required packages: open directory in terminal and type:

`pip install -r Requirements.txt`

For dash:
`pip install -r Requirements_dash.txt`

To install the package, open terminal in the directory and type:

`pip install -e .`      # . means this directory




Alternatively, the package can be installed directly from github, in this case the main folder is cloned to the local drive and named src
To install from github directly from Github:

Requirements:
`pip install -r https://raw.githubusercontent.com/PSenck/biomoni/main/Requirements.txt`

For dash: 
`pip install -r https://raw.githubusercontent.com/PSenck/biomoni/main/Requirements_dash.txt`

biomoni:
`pip install -e git+https://github.com/PSenck/biomoni.git#egg=biomoni`





