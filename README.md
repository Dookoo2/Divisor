# Private Key Recovery Tool

Funny high-performance C++ tool for recovering private keys from known public keys using a combination of distinguished points and parallel computation with divisors.
Like a joke, not for real work!

## Features

- **Distinguished Points Algorithm**: Efficient key recovery using trap-based searching
- **Multi-threaded**: Utilizes OpenMP for parallel computation
- **Bloom Filter**: Optimized memory usage with SIMD-accelerated bloom filters
- **Progress Monitoring**: Real-time statistics and progress reporting
- **Flexible Search**: Supports both forward and random search patterns

## Requirements

- C++17 compatible compiler
- OpenMP support
- x86_64 architecture with AVX2 support (for SIMD bloom filters)

## Build

```bash
g++ -O3 -march=native -fopenmp -std=c++17 -o keyrecover main.cpp


Sample 1
root@DESKTOP-BD9V01U:/mnt/e/Mark3# ./Mark3 -k 025e466e97ed0e7910d3d90ceb0332df48ddf67d456b9e7303b50a3d89de357336 -r 1:F02B35A358F -q 10000000 --order forward

=========== Phase-0: RAM summary ===========
DP table : 521MB  (10000000 / 13333334 slots, load 75.00%)
Bloom    : 19.1MB
Total    : 540MB

============= Phase-1: Working =============
Full time: 4.02 s
Speed    : 692.7kkey/s
Shrink   : 0
Offset   : 0x0000000000000000000000000000000000000000000000000000000000000000

============== Phase-2: RESULT =============
Divisor  : 0x00000000000000000000000000000000000000000000000000000000002a7a2d
Multiple : 0x000000000000000000000000000000000000000000000000000000000054702b
Offset   : 0x0000000000000000000000000000000000000000000000000000000000000000
Private  : 0x00000000000000000000000000000000000000000000000000000e02b35a358f
Full time: 4.02 s

Sample 2
root@DESKTOP-BD9V01U:/mnt/e/Mark2# ./Mark3 -k 03f46f41027bbf44fafd6b059091b900dad41e6845b2241dc3254c7cdd3c5a16c6 -r 1:42BD43C2E9354 -q 10000000 --order f
orward

=========== Phase-0: RAM summary ===========
DP table : 521MB  (10000000 / 13333334 slots, load 75.00%)
Bloom    : 19.1MB
Total    : 540MB

============= Phase-1: Working =============
Full time: 261.98 s
Speed    : 1.2Mkey/s
Shrink   : 4
Offset   : 0x00000000000000000000000000000000000000000000000000016f39f5dfc971

============== Phase-2: RESULT =============
Divisor  : 0x0000000000000000000000000000000000000000000000000000000001899953
Multiple : 0x00000000000000000000000000000000000000000000000000000000007aab31
Offset   : 0x00000000000000000000000000000000000000000000000000016f39f5dfc971
Private  : 0x00000000000000000000000000000000000000000000000000022bd43c2e9354
Full time: 261.98 s
 ```

## Getting Started
```bash
 g++ -std=c++17 -Ofast -funroll-loops -ftree-vectorize -fstrict-aliasing -fno-semantic-interposition -fvect-cost-model=unlimited -fno-trapping-math -fipa-ra -fipa-modref -flto -fassociative-math -fopenmp -mavx2 -mbmi2 -madx Mark3.cpp Int.cpp SECP256K1.cpp Point.cpp Random.cpp IntMod.cpp IntGroup.cpp Timer.cpp -o Mark3
 ```

## TIPS
BTC: bc1qtq4y9l9ajeyxq05ynq09z8p52xdmk4hqky9c8n
