{-# LANGUAGE ExplicitForAll #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}

module Partition
  ( greedyPart,
    getPartition,
    writePartition,
    getBlocksMatchGraph,
    longestSubstring,
  )
where

import Data.ByteString.Char8 qualified as BS
import Data.IntSet (IntSet)
import Data.IntSet qualified as IntSet
import Data.List qualified as List
import Data.Maybe (fromMaybe, mapMaybe)
import Genomes (GenesIRsF, GenesIRsR, Genome (..), Idx, IntergenicGenome (..), Matcher (..), RigidFlexibleReverseMatcher (..), decIdx, idxDist, idxToInt, incIdx, mkIdx, writeFGenome, writeIR, writeRGenome)
import LocalBase

type Breakpoint = [Idx]

data Partition g1 g2 where
  GenomePartition :: (Genome g1, Genome g2) => g1 -> g2 -> Breakpoint -> Breakpoint -> Partition g1 g2

instance Show (Partition GenesIRsR GenesIRsF) where
  show (GenomePartition g h bg bh) = unlines [combiStr subs_g, combiStr subs_h]
    where
      combiStr strs = unwords $ interleavelists strs (replicate (length strs - 1) "|")
      subs_g = zipWith (getSub show g) bg $ tail bg
      subs_h = zipWith (getSub show h) bh $ tail bh
      getSub write g' i succi = write $ subGenome (incIdx i) succi g'

getPartition :: GenesIRsR -> GenesIRsF -> Partition GenesIRsR GenesIRsF
getPartition = greedyPart RFRM

writePartition :: Partition GenesIRsR GenesIRsF -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString)
writePartition (GenomePartition g h bg bh) = (genes_bs_g, irs_bs_g, genes_bs_h, irs_bs_h)
  where
    genes_bs_g = combiBS $ map fst bssg
    irs_bs_g = combiBS $ interleavelists (map snd bssg) ir_breaks_g
    genes_bs_h = combiBS $ map fst bssh
    irs_bs_h = combiBS $ interleavelists (map snd bssh) ir_breaks_h
    combiBS bss = BS.unwords $ interleavelists bss (replicate (length bss - 1) "|")
    ir_breaks_g = map (writeIR . (`getIR` g)) (init . tail $ bg)
    ir_breaks_h = map (writeIR . (`getIR` h)) (init . tail $ bh)
    bssg = zipWith (getSub (writeRGenome False) g) bg $ tail bg
    bssh = zipWith (getSub (writeFGenome False) h) bh $ tail bh
    getSub write g' i succi = write $ subGenome (incIdx i) succi g'

getBlocksMatchGraph :: RigidFlexibleReverseMatcher GenesIRsR GenesIRsF -> Partition GenesIRsR GenesIRsF -> [[Int]]
getBlocksMatchGraph macher (GenomePartition g h bg bh) =
  do
    sub_g <- zipWith (getSub g) bg $ tail bg
    let sub_hs = zipWith (getSub h) bh $ tail bh
    return . map fst . filter (isMatch macher sub_g . snd) . zip [0 ..] $ sub_hs
  where
    getSub g' i succi = subGenome (incIdx i) succi g'

-- TODO: add withSin, see original code

-- | Greedy algorithm for string partition with intergenic regions or
-- flexible intergeic regions
greedyPart :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> g1 -> g2 -> Partition g1 g2
greedyPart matcher g h = GenomePartition g h (cleanList final_breaksG) (cleanList final_breaksH)
  where
    cleanList = map head . List.group . List.sort
    (final_breaksG, final_breaksH) = getAllLS IntSet.empty IntSet.empty [] []
    getAllLS genesOnBlocksG genesOnBlocksH breaksG breaksH =
      case longestSubstring matcher genesOnBlocksG genesOnBlocksH g h of
        Nothing -> (breaksG, breaksH)
        Just (((g_beg, g_end), (h_beg, h_end)), genesOnBlocksG', genesOnBlocksH') ->
          let breaksG' = (decIdx g_beg : g_end : breaksG)
              breaksH' = (decIdx h_beg : h_end : breaksH)
           in getAllLS genesOnBlocksG' genesOnBlocksH' breaksG' breaksH'

longestSubstring :: forall m g1 g2. (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> IntSet -> IntSet -> g1 -> g2 -> Maybe (((Idx, Idx), (Idx, Idx)), IntSet, IntSet)
-- Find longest substring of l1 and l2, ignore genes already in a block
-- returns the indices of the correspondent subgenomes in g and h (first and last gene)
longestSubstring matcher genesOnBlocksG genesOnBlocksH g h =
  case maybe_inds of
    Nothing -> Nothing
    Just inds@((g_beg, g_end), (h_beg, h_end)) ->
      Just
        ( inds,
          foldr (IntSet.insert . idxToInt) genesOnBlocksG [g_beg .. g_end],
          foldr (IntSet.insert . idxToInt) genesOnBlocksH [h_beg .. h_end]
        )
  where
    maybe_inds =
      maxWith (uncurry idxDist . fst)
        . mapMaybe (longestMatch . mkIdx)
        . filter (`IntSet.notMember` genesOnBlocksG)
        $ [1 .. size g]

    longestMatch :: Idx -> Maybe ((Idx, Idx), (Idx, Idx))
    -- find longest match between a subgenome of h and a suffix of g,
    -- ignore genes already in a block,
    -- returns the indices of the correspondent subgenomes in g and h (first and last gene)
    longestMatch ig = maxWith (uncurry idxDist . fst) $ do
      ih <- map mkIdx . filter (`IntSet.notMember` genesOnBlocksH) $ [1 .. size h]
      (endCPG, endCPH) <- fromMaybe [] . fmap (: []) $ commonPrefix ig ih
      return ((ig, endCPG), (ih, endCPH))

    commonPrefix :: Idx -> Idx -> Maybe (Idx, Idx)
    -- find common prefix of subgenome of g starting with ig
    -- and subgenome of h starting with ih
    -- returns Nothing if the prefix is empty
    -- returns indices for the ends of both subgenomes otherwise
    commonPrefix ig ih = if getGene ig g /= getGene ih h then Nothing else Just $ commonPrefix' ig ig ih ih

    commonPrefix' :: Idx -> Idx -> Idx -> Idx -> (Idx, Idx)
    commonPrefix' ig jg ih jh =
      if
          | jg == mkIdx (size g) -> (jg, jh)
          | jh == mkIdx (size h) -> (jg, jh)
          | idxToInt jg' `IntSet.member` genesOnBlocksG -> (jg, jh)
          | idxToInt jh' `IntSet.member` genesOnBlocksH -> (jg, jh)
          | not (isMatch matcher (subGenome jg jg' g) (subGenome jh jh' h)) -> (jg, jh)
          | otherwise -> commonPrefix' ig jg' ih jh'
      where
        jg' = incIdx jg
        jh' = incIdx jh