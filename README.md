![alt text](https://github.com/aleable/N-STARK-develop/blob/main/misc/logo.svg)

___

N-STARK (**N**on-**STA**tionary loads **R**outing on networ**K**s) is a Python implementation of the algorithms used in:

- [1] Alessandro Lonardi, Enrico Facca, Mario Putti and Caterina De Bacco. <i>Routing on networks with time dependent loads</i>. [arXiv] (ADD LINK).

This is a scheme capable of solving network routing problems with input loads injected in the nodes that change in time. Optimal solutions are computed solving a dynamical system of equations.

**If you use this code please cite [1].**

## What's included

- ```code```: contains the all the scripts necessary to run N-STARK, and a Jupyter notebook (```dashboard.ipynb```) with which it is possible to easily interact with them and visualize the results
- ```data/input```: contains the data needed to build the Bordeaux bus network, where the algorithm can be tested. The network topology has been pre-processed and extracted using [2]
- ```misc```: files used for the README.md
- ```setup.py```: setup file to build the Python environment

[2] Rainer Kujala, Christoffer Weckström, Richard K. Darst, Miloš N. Mladenović and Jari Saramäki, <a href="https://www.nature.com/articles/sdata201889">Scientific data <b>5</b>, 180089 (2018)</a>.<br/>

## How to use

To download this repository, copy and paste the following:

```bash
git clone https://github.com/aleable/N-STARK
```


**You are ready to test the code! But if you want to know how click [HERE](https://github.com/aleable/N-STARK/tree/main/code)**.

## Contacts

For any issues or questions, feel free to contact us sending an email to <a href="alessandro.lonardi@tuebingen.mpg.de">alessandro.lonardi@tuebingen.mpg.de</a>.

## License

Copyright (c) 2021 <a href="https://aleable.github.io/">Alessandro Lonardi</a> and <a href="https://www.cdebacco.com/">Caterina De Bacco</a>

<sub>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:</sub>

<sub>The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.</sub>

<sub>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.</sub>
