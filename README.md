**MINTsC** : Learning **M**ulti-way **INT**eractions from **s**ingle cell Hi-**C**
![ELECT diagram](/figures/intro.png)

## MINTsC Usage

### 1. Preparation


First of all, before cloning the MINTsC github package, go to the right directory that you would like to implement MINTsC. In the cmd terminal, do

```
cd {MINTsC directory}
```

then go on to the next step. {MINTsC directory} could be for example /storage10/kwangmoon/MINTsC.


#### 1. Repository clone

For cloning the github repository, again on the cmd terminal, run the linux code 

```
git clone git@github.com:keleslab/MINTsC.git
```

For R, we need the requirements as below : 


-   R: [R installation](https://www.r-project.org)  (>=4.2.2)

#### 2. Install/load required R packages

In R, run those codes that download the required packages for running MINTsC.

```
install.packages('pacman')
pacman::p_load(splines,dplyr,data.table,stringr,purrr,furrr,COUNT,optimParallel,
               RhpcBLASctl,purrr,doParallel,foreach)
```

Details about implementing codes can be found in the [Tutorial](https://github.com/keleslab/MINTsC/blob/main/code/scripts/Tutorial.ipynb) of this github, which uses the [Ramani et al.2017](https://www.nature.com/articles/nmeth.4155) scHi-C data as a tutorial.

