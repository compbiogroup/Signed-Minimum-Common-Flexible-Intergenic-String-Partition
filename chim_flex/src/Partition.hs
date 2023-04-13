{-# LANGUAGE ExplicitForAll #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}

module Partition
  ( greedyPart,
    getPartition,
    longestSubstring,
  )
where

import Data.IntSet (IntSet)
import Data.IntSet qualified as IntSet
import Data.List qualified as List
import Data.Maybe (fromMaybe, mapMaybe)
import Genomes (GenesIRsF, GenesIRsR, Genome (..), Idx, Matcher (..), RigidFlexibleMatcher (..), decIdx, idxDist, idxToInt, incIdx, mkIdx)
import LocalBase

type Breakpoint = [Idx]

data Partition g1 g2 where
  GenomePartition :: (Genome g1, Genome g2) => g1 -> g2 -> Breakpoint -> Breakpoint -> Partition g1 g2

getPartition :: GenesIRsR -> GenesIRsF -> Partition GenesIRsR GenesIRsF
getPartition = greedyPart RFM

reduced :: Partition GenesIRsR GenesIRsF -> (GenesIRsR, GenesIRsF)
reduced (Partition g h bg bh) = undefined

-- | Converts a partition into a pair of genomes representing a mapping into
-- a permutation compatible with the reduced genomes correspondent to the partition
mapToPerm :: Partition GenesIRsR GenesIRsF -> (GenesIRsR, GenesIRsF)
mapToPerm (Partition g h bg bh) = gr hr
  where
    (_, _, g_bal_ir) = toLists False g_bal
    (_, _, h_bal_ir) = toLists False h_bal
    gr = fromLists False sign gr_ls g_bal_ir
    hr = fromLists False sign hr_ls h_bal_ir
    (gr_ls, hr_ls) = (genomesToUniqueList ggs, genomesToUniqueList hhs)
    ggs = gseq part
    hhs = hseq part
    (sign, gmEmptyX) = case partType part of
      MCISP -> (Unsigned, gmEmpty)
      RMCISP -> (Signed, gmEmptyRev)

    -- Produce list of unique characters correspondent to genes
    genomesToUniqueList :: Seq Genome -> [Gene]
    genomesToUniqueList = genomesToUniqueList' [] m_orig . toList
    genomesToUniqueList' :: [Gene] -> GenomeMap GeneListForMap -> [Genome] -> [Gene]
    genomesToUniqueList' acc m [] = reverse acc
    genomesToUniqueList' acc m (g : gs) = genomesToUniqueList' ([v .. v + n - 1] ++ acc) m' gs
      where
        (_, lg, _) = toLists False g
        n = intToGene . coerce . genomeSize $ g
        v = head . unGeneListForMap . fromMaybe (error "Error on mapToPerm (1).") $ gmLookup g m
        m' =
          gmAlter
            g
            ( \old -> case old of
                Nothing -> error "Error on mapToPerm (2)."
                Just (MkGeneListForMap (x : xs)) -> MkGeneListForMap xs
            )
            m
    m_orig = fst $ foldl addGenome (gmEmptyX, 1 :: Gene) (hhs Seq.>< ggs)
    addGenome :: (GenomeMap GeneListForMap, Gene) -> Genome -> (GenomeMap GeneListForMap, Gene)
    addGenome (m, count) g = (m', count')
      where
        n = intToGene . coerce . genomeSize $ g
        count' = count + n
        m' =
          gmAlter
            g
            ( \old -> case old of
                Nothing -> MkGeneListForMap [count]
                Just (MkGeneListForMap l) -> MkGeneListForMap (count : l)
            )
            m

-- | Converts a partition into a pair of genomes representing the reduced genomes
-- correspondent to the partition
reduced :: Partition -> (Genome, Genome)
reduced (ValidPartition part) = (gr, hr)
  where
    gr = fromLists False sign gr_ls (toList $ gbps part)
    hr = fromLists False sign hr_ls (toList $ hbps part)
    (gr_ls, hr_ls) = (genomesToGenes ggs, genomesToGenes hhs)
    ggs = gseq part
    hhs = hseq part
    (sign, gmEmptyX) = case partType part of
      MCISP -> (Unsigned, gmEmpty)
      RMCISP -> (Signed, gmEmptyRev)

    genomesToGenes = map genomeToGene . toList
    genomeToGene g = fromMaybe (error "Error on reduced.") $ gmLookup g m
    m = fst $ foldl addGenome (gmEmptyX, 1 :: Gene) (hhs Seq.>< ggs)
    addGenome (m, count) g = (m', count')
      where
        count' = if keepOld then count else count + 1
        (m', _, keepOld) = gmLookupInsert g count m

-- TODO: add withSin, see original code
greedyPart :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> g1 -> g2 -> Partition g1 g2
greedyPart matcher g h = GenomePartition g h (List.sort final_breaksG) (List.sort final_breaksH)
  where
    (final_breaksG, final_breaksH) = getAllLS IntSet.empty IntSet.empty [] []
    getAllLS genesOnBlocksG genesOnBlocksH breaksG breaksH =
      case longestSubstring matcher genesOnBlocksG genesOnBlocksH g h of
        Nothing -> (breaksG, breaksH)
        Just (((g_beg, g_end), (h_beg, h_end)), genesOnBlocksG', genesOnBlocksH') ->
          let breaksG' = (decIdx g_beg : g_end : breaksG)
              breaksH' = (decIdx h_beg : h_end : breaksG)
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