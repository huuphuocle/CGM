# Class group method (CGM) for factorizing integers

## About
This project is started during my Master 2 Applied Algebra of University Paris-Saclay. It is an implementation of Class Group Method for integer factorization, i.e., given

`N = 6311189789556231825830216777410344573274773056312746064774869867452658507188`,

this program can write `N` into

`N = 2^2 * 3 * 43 * 1117 * 10055557348186662257983269247161283 * 1088935655927705848145865664713542363`.

## Installation
### Dependencies
This program only requires [GMP](https://gmplib.org/) to work. You can download the latest version of GMP (v 6.2.1) here [https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz](https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz).

### Build & Usage
To compile the program, you can use `makefile`: `make all`
To use the program:

* Create a file `<data-file>` whose the first line contains the number of primes need to be factorized and every following line contains a prime.
* Run `./cgm <data-file>`



