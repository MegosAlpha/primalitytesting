# Primality Testing
Companion to my Number Theory paper "A Brief History of Primality Testing", consisting of implementations of various primality tests in Python.

## Original Benchmarks
In the paper, I used my home PC for all benchmarking. It has a Ryzen 9 3900X and 32GB of RAM, so single-threaded performance was the most significant restraint.

## Implementations
- Baillie-PSW
- Agrawal-Kayal-Saxena (AKS): Requires a working installation of SageMath
- Miller-Rabin
- A pigeonhole principle extension to Miller-Rabin (not the Miller Test that relies on the Extended Riemann Hypothesis)
- Fermat Test
- Iterator version of the Sieve of Eratosthenes

## Requirements
Python 3.6 or greater, since I use f-strings. Also some modern version of SymPy, since I avoid implementing the Bernstein-Lenstra-Pila algorithm and the Jacobi Symbol. I used version 1.8 alongside Python 3.8.5 with Anaconda at time of writing. For AKS, you need SageMath -- I used version 9.2 on Python 3.8.8.
