# Installation

Linux only !

## Install through pip
```shell
$ pip install puckepy
```

## Local installation
Requires Python3 `>= 3.12` and Rust `>= 1.85`

Install **Rust**
```shell 
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Download the repository and enter the directory
```shell 
$ git clone https://github.com/jrihon/puckepy.git
$ cd puckepy/
```

Install the **maturin** framework
```shell
$ pip install maturin
```
Create a **virtual env** through pip 
```shell 
$ python3 -m venv .venv
$ source .venv/bin/activate
```

Compile the **puckepy library**
```shell
$ maturin develop
```

For other installation methods of maturin, visit the [maturin.rs website](https://www.maturin.rs/installation)
