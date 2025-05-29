# Applying COWBOYS to the benchmark of High-Dimensional Bayesian Optimization for discrete sequence optimization

[![Link to Project website](https://img.shields.io/badge/GitHub-Project_Website-100000?logo=github&logoColor=white)](https://machinelearninglifescience.github.io/hdbo_benchmark)
[![Link to Project website](https://img.shields.io/badge/GitHub-poli_docs-100000?logo=github&logoColor=white)](https://machinelearninglifescience.github.io/poli-docs)
[![Tests on hdbo (conda, python 3.10)](https://github.com/MachineLearningLifeScience/hdbo_benchmark/actions/workflows/tox-lint-and-pytest.yml/badge.svg)](https://github.com/MachineLearningLifeScience/hdbo_benchmark/actions/workflows/tox-lint-and-pytest.yml)

This repository contains the code to apply COWBOYS to a benchmark of **high-dimensional Bayesian optimization** over discrete sequences using [poli](https://github.com/MachineLearningLifeScience/poli) and [poli-baselines](https://github.com/MachineLearningLifeScience/poli-baselines).

### Running your solver locally

We provide a `requirements.txt`/`environment.yml` you can use to create an environment for running the benchmarks. Afterwards, install the additional packages:

```bash
conda create -n hdbo_benchmark python=3.10
conda activate hdbo_benchmark
pip install -r requirements.txt
pip install gauche
pip install selfies
pip install -e .
```
and change the WANDB_PROJECT and WANDB_ENTITY in src/hdbo_benchmark/utils/constants.py to link your WANDB account for easy visualisation of results.

To run COWBOYS across all the molecular search benchmark problems, simply run the bash script

```bash
./run.sh
```

assuming `hdbo_benchmark` is an environment in which you can run your solver, and in which this package is installed. The first time you run this script might take some time, while the testing environment is prepared.

Individual problems can be ran for specific seeds using run.py.
