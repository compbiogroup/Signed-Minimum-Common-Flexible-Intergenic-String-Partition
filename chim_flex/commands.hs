{-# LANGUAGE OverloadedStrings #-}

import LocalBase
import Genomes
import Partition
import PApprox

g = readRGenome True Signed "1 2" "1 1 2"
h = readRGenome True Signed "-2 -1" "0 1 3"