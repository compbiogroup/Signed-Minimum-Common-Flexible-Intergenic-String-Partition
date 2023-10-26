{-# LANGUAGE ExplicitForAll #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}

-- \|
-- Module      : PGreedy
-- Description : Implementation of the greedy algorithm for string partition adapted
-- to include intergenic region information.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module PGreedy
  ( greedyPartition,
    longestSubstring,
    commonPrefix,
  )
where

import Data.EnumSet (EnumSet)
import Data.EnumSet qualified as EnumSet
import Data.Foldable (toList)
import Data.Maybe (fromMaybe, mapMaybe, maybeToList)
import Genomes (Gene, Genome (alphabet, isGene, size, subGenome), Idx, Matcher (isDirectMatch, isReverseMatch), decIdx, idxDist, incIdx, mkIdx, positionMap, singletonOnBoth)
import LocalBase
import Partition (CommonPartition, mkCommonPartition2)

-- | Greedy algorithm for string partition with intergenic regions or
-- flexible intergeic regions. If withSingleton is True it looks first
-- for matched blocks that contain a singleton.
greedyPartition :: (Genome g1, Genome g2, Matcher m g1 g2) => Bool -> m g1 g2 -> g1 -> g2 -> CommonPartition g1 g2
greedyPartition withSingleton matcher g h = mkCommonPartition2 matcher g bg h bh
  where
    bg = EnumSet.delete (mkIdx $ size g) . EnumSet.delete 0 $ final_breaksG
    bh = EnumSet.delete (mkIdx $ size h) . EnumSet.delete 0 $ final_breaksH
    (final_breaksG, final_breaksH) = getAllLS all_singletons EnumSet.empty EnumSet.empty EnumSet.empty EnumSet.empty
    getAllLS singletons genesOnBlocksG genesOnBlocksH breaksG breaksH =
      case longestSubstring' singletons of
        (Nothing, Nothing) -> (breaksG, breaksH)
        (Just other_singletons, Nothing) -> getAllLS other_singletons genesOnBlocksG genesOnBlocksH breaksG breaksH
        (maybe_other_singletons, Just (((g_beg, g_end), (h_beg, h_end)), genesOnBlocksG', genesOnBlocksH')) ->
          let breaksG' = EnumSet.insert (decIdx g_beg) (EnumSet.insert g_end breaksG)
              breaksH' = EnumSet.insert (decIdx h_beg) (EnumSet.insert h_end breaksH)
           in getAllLS (fromMaybe [] maybe_other_singletons) genesOnBlocksG' genesOnBlocksH' breaksG' breaksH'
      where
        longestSubstring' (singleton : other_singletons) = (Just other_singletons, longestSubstring (Just singleton) matcher genesOnBlocksG genesOnBlocksH g h)
        longestSubstring' [] = (Nothing, longestSubstring Nothing matcher genesOnBlocksG genesOnBlocksH g h)

    all_singletons =
      if withSingleton
        then filter (singletonOnBoth posMapG posMapH) . toList $ alphabet g
        else []
    posMapG = positionMap g
    posMapH = positionMap h

longestSubstring :: forall m g1 g2. (Genome g1, Genome g2, Matcher m g1 g2) => Maybe Gene -> m g1 g2 -> EnumSet Idx -> EnumSet Idx -> g1 -> g2 -> Maybe (((Idx, Idx), (Idx, Idx)), EnumSet Idx, EnumSet Idx)
-- Find longest substring of l1 and l2, ignore genes already in a block
-- returns the indices of the correspondent subgenomes in g and h (first and last gene).
-- If singleton has a value (Just x), only consider blocks containing that gene x.
longestSubstring maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h =
  case maybe_inds of
    Nothing -> Nothing
    Just inds@((g_beg, g_end), (h_beg, h_end)) ->
      Just
        ( inds,
          foldr EnumSet.insert genesOnBlocksG [g_beg .. g_end],
          foldr EnumSet.insert genesOnBlocksH [h_beg .. h_end]
        )
  where
    maybe_inds =
      maxWith (uncurry idxDist . fst)
        . mapMaybe longestMatch
        . filter (`EnumSet.notMember` genesOnBlocksG)
        $ [1 .. mkIdx (size g)]

    longestMatch :: Idx -> Maybe ((Idx, Idx), (Idx, Idx))
    -- find longest match between a subgenome of h and a suffix of g,
    -- ignore genes already in a block,
    -- returns the indices of the correspondent subgenomes in g and h (first and last gene)
    longestMatch ig = maxWith (uncurry idxDist . fst) $ do
      ih <- filter (`EnumSet.notMember` genesOnBlocksH) $ [1 .. mkIdx (size h)]
      maybeToList $ commonPrefix maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h ig ih

commonPrefix :: (Matcher m g1 g2, Genome g1, Genome g2) => Maybe Gene -> m g1 g2 -> EnumSet Idx -> EnumSet Idx -> g1 -> g2 -> Idx -> Idx -> Maybe ((Idx, Idx), (Idx, Idx))
-- find common prefix of subgenome of g starting with ig
-- and subgenome of h starting with ih
-- or if the elements is ig and ih are a reversed match
-- find common prefix of subgenome of g starting with ig
-- and subgenome of rev(h) starting with ih
-- returns Nothing if the prefix is empty
-- returns indices for the ends of both subgenomes otherwise
commonPrefix maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h ig ih = do
  (endCPG, endCPH, rev) <-
    if
        | isDirectMatch matcher (subGenome ig ig g) (subGenome ih ih h) -> Just $ commonPrefix' False ig ih
        | isReverseMatch matcher (subGenome ig ig g) (subGenome ih ih h) -> Just $ commonPrefix' True ig ih
        | otherwise -> Nothing
  case maybe_singleton of
    Just singleton ->
      if singleton `isGene` subGenome ig endCPG g
        then return ((ig, endCPG), if rev then (endCPH, ih) else (ih, endCPH))
        else Nothing
    Nothing -> return ((ig, endCPG), if rev then (endCPH, ih) else (ih, endCPH))
  where
    commonPrefix' :: Bool -> Idx -> Idx -> (Idx, Idx, Bool)
    commonPrefix' rev jg jh =
      if
          | jg == mkIdx (size g) -> (jg, jh, rev)
          | jh_end > mkIdx (size h) -> (jg, jh, rev)
          | jh_beg < 0 -> (jg, jh, rev)
          | jg' `EnumSet.member` genesOnBlocksG -> (jg, jh, rev)
          | jh' `EnumSet.member` genesOnBlocksH -> (jg, jh, rev)
          | not (testMatch matcher (subGenome jg jg' g) (subGenome jh_beg jh_end h)) -> (jg, jh, rev)
          | otherwise -> commonPrefix' rev jg' jh'
      where
        jg' = incIdx jg
        (jh_beg, jh_end, jh') = if rev then (decIdx jh, jh, decIdx jh) else (jh, incIdx jh, incIdx jh)
        testMatch = if rev then isReverseMatch else isDirectMatch
