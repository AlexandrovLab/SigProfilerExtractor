
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.6] - 2026-01-06

### Changed
- Updated default COSMIC version from 3.4 to 3.5. Added support for COSMIC v3.5 signatures in the `cosmic_version` parameter.
- Updated SigProfilerAssignment dependency requirement from >=1.0.1 to >=1.1.0 to support COSMIC v3.5 signatures.

## [1.2.5] - 2025-10-28

### Added
- Implemented a CI/CD pipeline with Travis CI to automate the building and publishing of Docker images to Docker Hub.
- Added a Dockerfile to the repository for containerization. Documentation on how to use the Dockerfile needs to be added to the README.

## [1.2.4] - 2025-10-20

### Added
- Added the `assignment_cpu` parameter to independently control the number of CPU cores used for the signature assignment step. This change enables full support for the parallel processing enhancements in **SigProfilerAssignment v1.0.0**, allowing for significant performance improvements and more granular resource control.

## [1.2.3] - 2025-09-19

### Added
- Added support for rn7 and mm39 genomes in SigProfilerExtractor.

## [1.2.2] - 2025-08-11

### Added
- Added mutation count and stability to 4608 plots in All_Solutions
- Add a stop parameter to the CLI to stop after de novo extraction.

## [1.2.1] - 2025-05-14

### Fixed
- Fixed an issue where the CLI was returning a non-zero exit code when the `--help` flag was passed.

 ### Added
- Added `pyproject.toml` for modern Python packaging support.

## [1.2.0] - 2025-02-11

### Changed
- Updated dependencies: Now requires **Pandas >= 2.0.0**, **NumPy >= 2.0.0**, and **Python >= 3.9**.
- Dropped support for **Python 3.8**
- **Intel-based MacBooks are no longer supported** due to upstream changes in **PyTorch**, which has dropped support for macOS x86_64. Users with Intel-based MacBooks will need to migrate to Apple Silicon (M1/M2) or use a Linux-based development environment.

## [1.1.25] - 2024-12-09

### Added
- Introduced a Command-Line Interface (CLI) for SigProfilerExtractor, enabling users to interact with the tool via terminal commands.

### Updated
- Improved the formatting of the parameter table for sigProfilerExtractor function for better readability and consistency.
- The CI/CD badge link has been fixed.
