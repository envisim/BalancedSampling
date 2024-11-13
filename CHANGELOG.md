# Changelog

All notable changes to the BalancedSampling package will be documented in this
file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed
- fixed bug in `cube`.
- fixed bug in `IndexList`.

## [2.0.6] - 2024-02-15
### Fixed
- `cube, lcube, cubestratified, lcubestratified` error regarding calculation of vectors in null space.
- `cubestratified, lcubestratified` now slightly more time-efficient for large data sets.

## [2.0.5] - 2024-02-09
### Added
- Examples for `lpm1, rpm, spm, lcps`.

### Fixed
- Bugfixes (`vsb`)

## [2.0.4] - 2024-02-07
### Fixed
- `cubestratified`, `lcubestratified`, breaking error
- `lcube`, bug impacting cases when initial probabilities were <= 0 or >= 1

## [2.0.3] - 2024-02-06
### Fixed
- Bugfixes (`lcubestratified`)

## [2.0.2] - 2024-02-02
- Removed SystemRequirements
- Added .Rbuildignore

## [2.0.1] - 2024-01-19
- Updated Rcpp dependency to 1.0.12.

## [2.0.0] - 2023-11-24
All functions have been rewritten, thus seeds are not respected between this version and previous versions.

### Added
- `lpm1s`, an efficient variant of `lpm1`, equivalent for populations without any equally distant neighbours.
- `genpopUniform` generates a uniform population.
- `genpopPoisson` generates a poisson cluster process.
- `sblb` gives a spatial balance measure using a local balance approach.

### Changed
- Parameter `Xbal` changed to `x` for `cube, cubestratified`.
- Parameter `eps` added for comparison between floats for `cube, cubestratified, hlpm2, lpm1, lpm2, spm, rpm, scps, lcps, lcube, lcubestratified`.
- `hlpm` now `hlpm2`. Parameter `p` changed to `prob`, `X` changed to `x`.
- `hlpm2, lpm1, lpm2, lcube, lcubestratified, scps, lcps, sb` uses a kd-tree to find nearest neighbour. These methods now accepts the parameter `type`, setting the splitting method of the kd-tree, and the parameter `bucketSize`, setting the maximal size of the buckets of the kd-tree.
- `lpm` is now an alias for `lpm2`.
- `probabilities` changed to `getPips`. Parameter `a` changed to `x`.
- Parameter `s` changed to `sample` for `sb`.
- Parameter `rand` added for `scps`, allowing for permanent random numbers.
- The algorithm for `scps` now chooses a unit at random, intead of traversing the list in order. If parameter `rand` is used, a fixed order is applied.
- Parameter `k` added for `vsb`, deciding on number of neighbours to construct means around.

### Removed
- `flightphase, landingphase, lcubeflightphase, lcubelandingphase`.
- Local Pivotal Method `lpm` using a window search.
- `scps_coord` removed, now accessed through optimal parameter in `scps`.
- `scps_getrand` removed.

## [1.6.3] - 2022-06-29

### Added
- Locally correlated Poisson sampling method `lcps`.

## [1.5.4] - 2018-09-03

### Added
- Added hierarchical local pivotal method.

### Removed
- lpm2_kdtree moved to new package SamplingBigData.

## [1.5.1] - 2016-01-27

### Added
- Added a faster kd-tree implementation of lpm2 by Jonathan Lisic.

## [1.4.0] - 2014-04-14

### Fixed
- Fixed a numerical stability issue in the cube-based methods.

## [1.3.0] - 2014-04-07

### Fixed
- Fixed a serious bug in cubestratified and lcubestratified.

## [1.2.0] - 2014-03-31

### Added
- Added stratified balanced sampling (cubestratified) and stratified doubly balanced sampling (lcubestratified).

## [1.1.0] - 2014-03-26

### Added
- Added separate functions for the flight phase and the landing phase of the cube method and the local cube method.

### Changed
- Renamed inclusionprobabilities to probabilities to not have a conflict with the package "sampling".
- Improved the speed of scps.

## [1.0.0] - 2014-03-19
- First release of version 1.0
