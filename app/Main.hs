{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}

-- |
-- Module      : MCFISP2K
-- Description : Algorithms for Genome partition problems with flexible intergenic regions.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module Main (main) where

import Control.Concurrent.ParallelIO.Global (parallel_, stopGlobalPool)
import Control.DeepSeq (force)
import Data.ByteString.Builder (intDec, toLazyByteString)
import Data.ByteString.Char8 qualified as BS
import Data.ByteString.Lazy qualified as LBS
import Data.Time (diffUTCTime, getCurrentTime)
import Genomes (GenesIRsF, GenesIRsR, RigidFlexibleReverseMatcher (..), ChromType(Linear), Sign (..), readFGenome, readRGenome)
import LocalBase
import Options.Applicative
import PApprox (approxPartition, approxLowerBound)
import PFpt (fptPartition)
import PGreedy (greedyPartition)
import PSOAR (soarPartition)
import Partition (CommonPartition, getBlocksCorrespondence, writePartition)
import Text.Printf (printf)

-- TODO: Make necessary adaptations to include multiple and circular chromosomes in the partition algorithms

data Args = Args
  { partAlg :: PartitionAlgorithm,
    input :: String,
    output :: String,
    noParallel :: Bool,
    signed :: Sign
  }

data PartitionAlgorithm = SOAR | SOARComb | Greedy | GreedySin | Approx | FPT | ApproxLB deriving (Read, Show)

argsParser :: Parser Args
argsParser =
  Args
    <$> option
      auto
      ( long "algorithm"
          <> short 'a'
          <> metavar "ALGO"
          <> help "Algorithm to use for the partition problem. The options are SOAR, Greedy, GreedySin, Approx, and FPT. There is also and ApproxLB option, which produces and integer corresponding to a lower bound of the number of breakpoints in the common partition (counting both partitions) instead of a partition on the output."
      )
    <*> strOption
      ( long "input"
          <> short 'i'
          <> metavar "IFILE"
          <> help "Input file. Each 4 lines of the input file correspond to a instance, each line has a list of comma or space separated values, and represent in order the origin string, the origin intergenic region list, the target string, and the target intergenic region list."
      )
    <*> strOption
      ( long "outfile"
          <> short 'o'
          <> metavar "OFILE"
          <> help "Output file. For each instance five lines are produces in the file, the first four lines correspond to two reduced genomes produced with the partition algorithm (each character of the string correspond to a block of the partition and each integer of the intergenic regions list correspond to a breakpoint). The last line shows the wall clock time required to produce the partition."
      )
    <*> switch
      ( long "no-par"
          <> help "Do not process the genomes in parallel."
      )
    <*> flag
      Unsigned
      Signed
      ( long "signed"
          <> short 's'
          <> help "Whether the input Strings are signed."
      )

opts :: ParserInfo Args
opts =
  info
    (argsParser <**> helper)
    ( fullDesc
        <> progDesc
          "Algorithm for genome partition problems with flexible intergenic regions."
    )

main :: IO ()
main = do
  args <- execParser opts
  contents <- BS.readFile (input args)
  let quadruples = zip [(1 :: Int) ..] . toQuadruples . filter ((/= '#') . BS.head) . BS.lines $ contents
  if noParallel args
    then mapM_ (runOne args) quadruples
    else do
      parallel_ $ map (runOne args) quadruples
      stopGlobalPool
  where
    runOne args (i, bstrs) = do
      start <- getCurrentTime
      !bstrs' <- fmap force (produceBlockMatch (partAlg args) (signed args) bstrs)
      end <- getCurrentTime
      let time = BS.pack . (show :: Double -> String) . realToFrac $ diffUTCTime end start
      BS.writeFile (output args ++ "_" ++ printf "%04d" i) . BS.unlines $ bstrs' ++ ["# Time: " <> (BS.pack . show $ time)]

    toQuadruples (s1 : i1 : s2 : i2 : ss) = (s1, i1, s2, i2) : toQuadruples ss
    toQuadruples [] = []
    toQuadruples _ = error "Incorrect number of lines."

produceBlockMatch :: PartitionAlgorithm -> Sign -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString) -> IO [BS.ByteString]
produceBlockMatch alg sign (s1, i1, s2, i2) =
  case alg of
    ApproxLB -> return [BS.pack . show $ approxLowerBound RFRM g h]
    _ -> do
      (maybe_part, fullComp) <- getPartition alg g h
      case maybe_part of
        Nothing -> if fullComp then return ["# Error: Execution finish but no partition was found"] else return ["# No solution found in time"]
        Just part ->
          let (s1', i1', s2', i2') = writePartition part
              bc = writeBlocksCorrespondence $ getBlocksCorrespondence RFRM part
           in return [s1', i1', s2', i2', bc, if fullComp then "# Exact solution" else "# Partial solution"]
  where
    g = readRGenome Linear True sign s1 i1
    h = readFGenome Linear True sign s2 i2

writeBlocksCorrespondence :: [[Int]] -> BS.ByteString
writeBlocksCorrespondence = BS.unwords . (\l -> interleavelists l (replicate (length l - 1) "|")) . map (BS.unwords . map (LBS.toStrict . toLazyByteString . intDec))

getPartition :: PartitionAlgorithm -> GenesIRsR -> GenesIRsF -> IO (Maybe (CommonPartition GenesIRsR GenesIRsF), Bool)
getPartition FPT g h = fptPartition 600000000 RFRM g h
getPartition SOAR g h = return . (,True) . Just $ soarPartition False RFRM g h
getPartition SOARComb g h = return . (,True) . Just $ soarPartition True RFRM g h
getPartition Greedy g h = return . (,True) . Just $ greedyPartition False RFRM g h
getPartition GreedySin g h = return . (,True) . Just $ greedyPartition True RFRM g h
getPartition Approx g h = return . (,True) . Just $ approxPartition RFRM g h
getPartition _ _ _ = error patternError
