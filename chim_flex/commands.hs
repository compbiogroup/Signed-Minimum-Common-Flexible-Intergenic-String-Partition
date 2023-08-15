{-# LANGUAGE OverloadedStrings #-}

import LocalBase
import Genomes
import Partition
import PApprox
import PFpt
import PSOAR

g = readRGenome True Signed "-1 -39 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1" "0 0 0 0 0 0 0 0 0 0 0 0 0 0"
h = readRGenome True Signed "-1  -1  1 39  1 -1 -1 -1 -1 -1 -1 -1 -1" "0 0 0 0 0 0 0 0 0 0 0 0 0 0"