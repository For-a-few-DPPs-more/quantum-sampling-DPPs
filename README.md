# quantum-sampling-DPPs
We provide here notebooks associated with the manuscript 

[On sampling determinantal and Pfaffian point processes on a quantum computer]()

by [Rémi Bardenet](https://rbardenet.github.io/), [Michaël Fanuel](https://mrfanuel.github.io/) and Alexandre Feller.


## Dependencies
Python >= 3.8.2

For sampling Determinantal PPs, please install
- [qiskit](https://qiskit.org/) 
- [numpy](https://numpy.org/) 
- [matplotlib](https://matplotlib.org/)
- [pandas](https://pandas.pydata.org/)

For sampling Pfaffian PPs, also install
- [pfapack](https://pypi.org/project/pfapack/)
- qiskit-nature: please use the fork of qiskit-nature available at
[`https://github.com/mrfanuel/qiskit-nature`](https://github.com/mrfanuel/qiskit-nature)


## Jupyter notebooks for reproducing the paper figures

- [slater_determinants.ipynb](https://github.com/For-a-few-DPPs-more/quantum-sampling-DPPs/blob/main/notebooks/slater_determinants.ipynb)
reproduces the simulations of projection DPPs.

- [pfaffian_sampling.ipynb](https://github.com/For-a-few-DPPs-more/quantum-sampling-DPPs/blob/main/notebooks/pfaffian_sampling.ipynb)
reproduces simulations of  Pfaffian  PPs.

