# sts-online experiments

This simulation pipeline was used in [Effective Online Bayesian Phylogenetics Via Sequential Monte Carlo With Guided Proposals](https://doi.org/10.1101/145219) using the following programs.

| Program    | Version |
|----------- | --------|
| sts-online | [v0.1](https://github.com/OnlinePhylo/sts/releases/tag/v0.1)           |
| MrBayes    | [v3.2.6](https://sourceforge.net/projects/mrbayes/files/mrbayes/3.2.6/)|
| PhyML      | [v3.2](https://github.com/stephaneguindon/phyml/releases/tag/v3.2.0)   |
| BAli-Phy   | [v2.3.8](https://github.com/bredelings/BAli-Phy/releases/tag/2.3.8)    |

# Requirements

## [sts-online](https://github.com/OnlinePhylo/sts)

The pipeline expects `sts-online` to be in the PATH.

## [MrBayes](http://mrbayes.sourceforge.net)

The pipeline expects `mb` to be in the PATH.

## [PhyML](https://github.com/stephaneguindon/phyml)

The pipeline expects `phyml` to be in the PATH.

## [BAli-Phy](http://bali-phy.org) package 

The pipeline requires the `trees-bootstrap` program to generate ASDSFs.
The pipeline expects `trees-bootstrap` to be in the PATH.

## virtualenv

``` shell
# install virtualenv if not already installed
pip install virtualenv

# create a virtualenv
virtualenv venv

# activate virtualenv
source venv/bin/activate

# install packages
pip install -r requirements.txt
```


## R packages

``` shell
# restore packrat project
R --no-restore --slave -e "0" --args --bootstrap-packrat
```

packrat installs the R packages into the `packrat/lib` directory.


## Bio++ and BppSuite

Running the simulations requires the tools `bppseqgen` from BppSuite 2.2.0 is required.

On Debian/Ubuntu:

``` shell
sudo apt-get install bppsuite
```

If BppSuite isn't available through your package manager, it can be installed from [source](http://biopp.univ-montp2.fr/repos/sources/bppsuite/bppsuite-2.2.0.tar.gz) into the `venv` prefix, so it'll be in your `PATH` when the virtualenv is activated:

``` shell
wget http://biopp.univ-montp2.fr/repos/sources/bppsuite/bppsuite-2.2.0.tar.gz
tar xf bppsuite-2.2.0.tar.gz
cd bppsuite-2.2.0
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../venv
make && make install
cd ..
```

Running the simulations requires the `trees-bootstrap` program from the BAli-Phy package

# Running the simulations

``` shell
cd comparison_to_mrbayes
python run_simulations.py
```

The simulations will probably take several days to complete. Multiple csv files will be produced in the `output` folder.


# Parsing results

```shell
cd comparison_to_mrbayes
Rscript --no-restore --slave -e rmarkdown::render("ess.Rmd")
```

The script will generate the file `sts.pdf`.
