# Longest-Common-Subsequence

Longest Common Subsequence dynamic programming algorithm.  

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
- [Acknowledgments](#acknowledgments)

## Introduction

This program efficiently finds the longest common subsequence of characters within two different arbitrarily large strings.
These algorithms are implemented based off of the Needleman Wunsch global sequence alignment and the Smith Waterman local
sequence alignment.

## Features

The algorithms both use dynamic programming to significantly decrease the run time. To brute force this alforithm the runtime would be O(2^m * n) assuming m <= n.
Both local and global sequence alignments have a runtime of O(n*m).

## Getting Started

To run the program you need to download three files: LongestCommonSubsequence.py, eyeless_genome_sequence.fasta, and PAX6_genomic_sequence.fast.
Store the files in the same directory and run them using the bash command: $ python3 ./LongestCommonSubsequence.py

### Prerequisites

This program is in python. So python is required as a well as numpy to run this program.

## Acknowledgments

My Professor Jessen Havill and peer Jacob Piskadlo assisted in implementing the alignment algorithms and back-tracking algorithms.

