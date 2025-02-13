
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
