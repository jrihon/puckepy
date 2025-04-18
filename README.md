# Puckepy

A library to describe and characterise molecules.

- `puckepy.formalism` to describe molecules quantitatively through various puckering formalisms
- `puckepy.confsampling` to provide the functionality from [pucke.rs](https://github.com/jrihon/puckers) in `Python`
- `puckepy.geometry` to describe molecules by elementary geometrical attributes

## Documentation
Consult the API documentation [here](https://github.com/jrihon/puckepy/blob/main/docs/documentation.md) !

<!--- The Python library has been annotated with python stub files `.pyi` to help your `LSP` with `definition on hover` functionality. -->

## Installation
+ Available for linux on `python3 --version >=3.12`. Consider making a new `.venv` or `conda` environment if your python version is older.
+ Requires no python dependencies, so it is safe to install in any `base` environment.

```shell
# pip venv
$ python3.12 -m venv .venv
$ source activate .venv
```

```shell
# Conda
$ conda create --name myenv python=3.12
$ conda activate myenv
```

Install `puckepy` in the environment.
```shell
$ pip install puckepy
```

If you are not on Linux, you might benefit from a local installation method. If you are on Windows, I suggest using Windows Subsystem for Linux (WSL2). \
Local installation protocol can be found [here](https://github.com/jrihon/puckepy/blob/main/docs/installation.md) !

## Citation 
If you've used pucke.py or any of its parts, please cite the following published article : 

[1] Rihon, J., Reynders, S., Bernardes Pinheiro, V. et al. The pucke.rs toolkit to facilitate sampling the conformational space of biomolecular monomers. J Cheminform 17, 53 (2025). https://doi.org/10.1186/s13321-025-00977-7

```bibtex
@article{Rihon_2025,
    title={The pucke.rs toolkit to facilitate sampling the conformational space of biomolecular monomers},
    volume={17},
    ISSN={1758-2946},
    url={http://dx.doi.org/10.1186/s13321-025-00977-7},
    DOI={10.1186/s13321-025-00977-7},
    number={1},
    journal={Journal of Cheminformatics},
    publisher={Springer Science and Business Media LLC},
    author={Rihon, Jérôme and Reynders, Sten and Bernardes Pinheiro, Vitor and Lescrinier, Eveline},
    year={2025},
    month=apr
}
```

## Author
Written and designed by [Jérôme Rihon](https://github.com/jrihon/jrihon)
