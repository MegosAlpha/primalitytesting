# Example benchmark
# This runs the first and second benchmarks from the paper.
# Notable exclusions are SymPy (since it was tested seperately) and
# AKS (since it requires SageMath). Note that the correct number of primes for the
# test is 1229, as Baillie-PSW, Sieve, and Pigeonhole Miller-Rabin ought
# to agree upon.

from primalitytesting import *
import timeit

start_time = timeit.default_timer()
# Run method here
k = SieveOfEratosthenes()
# Last prime less than 10000
while (k.next() < 9973):
    pass
print(f"Sieve: Runtime {timeit.default_timer()-start_time}s, Primes: {len(k.primes)}")

start_time = timeit.default_timer()
# k = 1
ftrp = set([n if fermat_test_random(n,1) else 2 for n in range(2,10000)])
print(f"Random k=1 Fermat Test: Runtime {timeit.default_timer()-start_time}s, Primes: {len(ftrp)}")

start_time = timeit.default_timer()
# k = 4
pmr = set([n if miller_rabin(n, 4) else 2 for n in range(2,10000)])
print(f"Miller-Rabin k=4: Runtime {timeit.default_timer()-start_time}s, Primes: {len(pmr)}")

start_time = timeit.default_timer()
pbpsw = set([n if baillie_psw(n) else 2 for n in range(2,10000)])
print(f"Baillie-PSW: Runtime {timeit.default_timer()-start_time}s, Primes: {len(pbpsw)}")

# Warning! This one is a bit longer.
start_time = timeit.default_timer()
pmrd = set([n if miller_rabin_deterministic(n) else 2 for n in range(2,10000)])
print(f"Pigeonhole Miller-Rabin: Runtime {timeit.default_timer()-start_time}s, Primes: {len(pmrd)}")

start_time = timeit.default_timer()
# k = 1
ftrp = fermat_test_random(10007,1)
print(timeit.default_timer()-start_time)

start_time = timeit.default_timer()
# k = 4
pmr = miller_rabin(10007, 4)
print(timeit.default_timer()-start_time)

start_time = timeit.default_timer()
pmrd = miller_rabin_deterministic(10007)
print(timeit.default_timer()-start_time)

start_time = timeit.default_timer()
# Disable cheating for Baillie-PSW by wiping the cache
lsU_cache = {}
pbpsw = baillie_psw(10007)
print(timeit.default_timer()-start_time)
