{-# LANGUAGE OverloadedStrings #-}

import LocalBase
import Genomes
import Partition

g = readRGenome True Signed "-1 -1 -1 -12 -1 -8 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1" "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
h = readRGenome True Signed "1 1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 -1 -12 2 -8 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1" "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"