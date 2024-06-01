#!/usr/bin/env bash

# get number of duplex reads in a fastq file passed as positional arg

faster --regex_string ".+;.+" "$@" | grep '^@' | wc -l