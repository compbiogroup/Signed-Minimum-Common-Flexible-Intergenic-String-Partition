{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE OverloadedStrings #-}

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
import Genomes (RigidFlexibleReverseMatcher (..), Sign (..), readFGenome, readRGenome, GenesIRsR, GenesIRsF)
import LocalBase
import Options.Applicative
import Partition (getBlocksCorrespondence, writePartition, CommonPartition)
import Text.Printf (printf)
import PFpt (fptPartition)

data Args = Args
  { input :: String,
    output :: String,
    noParallel :: Bool,
    signed :: Sign
  }

argsParser :: Parser Args
argsParser =
  Args
    <$> strOption
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
      !bstrs' <- produceBlockMatch (signed args) bstrs
      end <- getCurrentTime
      let time = BS.pack . (show :: Double -> String) . realToFrac $ diffUTCTime end start
      BS.writeFile (output args ++ "_" ++ printf "%04d" i) . BS.unlines $ bstrs' ++ ["# Time: " <> (BS.pack . show $ time)]

    toQuadruples (s1 : i1 : s2 : i2 : ss) = (s1, i1, s2, i2) : toQuadruples ss
    toQuadruples [] = []
    toQuadruples _ = error "Incorrect number of lines."

produceBlockMatch :: Sign -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString) -> IO [BS.ByteString]
produceBlockMatch sign (s1, i1, s2, i2) = do
  maybe_part <- getPartition g h
  case maybe_part of
    Nothing -> return ["# No solution found in time"]
    Just (part, fullComp) ->
      let (s1', i1', s2', i2') = writePartition part
          bc = writeBlocksCorrespondence $ getBlocksCorrespondence RFRM part
      in return [s1', i1', s2', i2', bc, if fullComp then "# Exact solution" else "# Partial solution"]
  where
    g = readRGenome True sign s1 i1
    h = readFGenome True sign s2 i2

writeBlocksCorrespondence :: [[Int]] -> BS.ByteString
writeBlocksCorrespondence = BS.unwords . (\l -> interleavelists l (replicate (length l - 1) "|")) . map (BS.unwords . map (LBS.toStrict . toLazyByteString . intDec))

getPartition :: GenesIRsR -> GenesIRsF -> IO (Maybe (CommonPartition GenesIRsR GenesIRsF, Bool))
getPartition = fptPartition 600000 RFRM
