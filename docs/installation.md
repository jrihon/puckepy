# Installation

## Pip install (Work In Progress)

Will be published to pip once the article is out


## Local installation
For other installation methods, visit the [maturin.rs website](https://www.maturin.rs/installation)

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
