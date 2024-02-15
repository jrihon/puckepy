Installation
==============


Conda installation
------------------

Local installation
------------------
Install **Rust**
::

    $ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Install **conda**
Go to conda_ website.
.. _conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html


Install the **maturin** framework
::

    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict
    $ conda install maturin


Create and activate a new conda environment
::

    $ conda create -n puckepy
    $ conda activate puckepy

Compile the puckepy library
::

    $ maturin develop

Every time you want to use the **puckepy** module, just call  
::

    $ conda activate puckepy
