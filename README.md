
# Buckyball Challenge

We use `buckyballgenerator.py` to generator the buckyball, including edges, nodes and neighbors.
## Algorithms:

### 1. Tensor network
`contraction2.m`
`tensor_product.m` is the auxiliary function computing the tensor contraction.
Tensor network gives the most accurate solution.
### 2. Monte Carlo simulation
`MCMC.m` and ` montecarloalgorithm.py`.
We exhaustively search for all possible ground state configurations. After about 1 million times sampling we obtain configuration number.

We use similar exhaustive sampling to compute partition function Z.
### 3. Wang-Landau Algorithm

`wanglandau.py` and `Wang.m`.
