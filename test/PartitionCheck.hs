{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE TemplateHaskell #-}

module PartitionCheck (tests) where

import Data.IntMap qualified as IntMap
import Data.EnumSet qualified as EnumSet
import Data.List qualified as List
import Data.Maybe (fromJust)
import Control.Arrow (first)
import Genomes
import GenomesCheck (GenomeWrapper (..), genGenome, genRGenome, rearrangeAndFlexibilizeGenome)
import Hedgehog
import Hedgehog.Gen qualified as Gen
import Hedgehog.Range qualified as Range
import Partition
import PSOAR
import PGreedy
import PApprox
import PFpt

prop_partitionConstructionMethodsAreEquivalent :: Property
prop_partitionConstructionMethodsAreEquivalent =
  property $ do
    g <- forAll (genRGenome 100)
    bps <- forAll (Gen.subsequence [1 .. mkIdx (size g - 1)])
    let bps' = 0 : bps ++ [mkIdx (Genomes.size g)]
        bs = zipWith (\i succi -> subGenome (incIdx i) succi g) bps' $ tail bps'
    assert $ mkPartitionFromBreakpoints g (EnumSet.fromList bps) == mkPartitionFromBlocks g bs

prop_commonPrefixFromBothAreEqual :: Property
prop_commonPrefixFromBothAreEqual =
  property $ do
    g <- forAll (genRGenome 100)
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h <- forAll $ rearrangeAndFlexibilizeGenome k g
    idx_g <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
    idx_h <- forAll . fmap head . Gen.shuffle . fromJust $ geneMapLookup (getGene idx_g g) (positionMap h)
    let ((g_beg, g_end), (h_beg, h_end)) = fromJust $ commonPrefix Nothing RFRM EnumSet.empty EnumSet.empty g h idx_g idx_h
    subG <- forAll . return $ subGenome g_beg g_end g
    subH <- forAll . return $ subGenome h_beg h_end h
    assert $ isMatch RFRM subG subH

prop_longestSubstringFromBothAreEqual :: Property
prop_longestSubstringFromBothAreEqual =
  property $ do
    g <- forAll (genRGenome 100)
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h <- forAll $ rearrangeAndFlexibilizeGenome k g
    let (((g_beg, g_end), (h_beg, h_end)), _, _) = fromJust $ longestSubstring Nothing RFRM EnumSet.empty EnumSet.empty g h
    let subG = subGenome g_beg g_end g
        subH = subGenome h_beg h_end h
    assert $ isMatch RFRM subG subH

prop_suboptimalRuleIntervalTurnsReplicasIntoSingletons :: Property
prop_suboptimalRuleIntervalTurnsReplicasIntoSingletons = suboptimalRuleTurnsReplicasIntoSingletons suboptimalRuleInterval

prop_suboptimalRulePairsTurnsReplicasIntoSingletons :: Property
prop_suboptimalRulePairsTurnsReplicasIntoSingletons = suboptimalRuleTurnsReplicasIntoSingletons suboptimalRulePairs

suboptimalRuleTurnsReplicasIntoSingletons :: (RigidFlexibleReverseMatcher GenesIRsR GenesIRsF -> GenesIRsR -> GenesIRsF -> (GenesIRsR, GenesIRsF)) -> Property
suboptimalRuleTurnsReplicasIntoSingletons rule =
  property $ do
    g <- forAll (genRGenome 100)
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h <- forAll $ rearrangeAndFlexibilizeGenome k g
    (g', h') <- forAll . return $ rule RFRM g h
    let sigG = getSingletons g
        sigH = getSingletons h
        sigG' = getSingletons g'
        sigH' = getSingletons h'
    assert $ sigG `List.isSubsequenceOf` sigG'
    assert $ sigH `List.isSubsequenceOf` sigH'
  where
    getSingletons k = List.sort $ foldr (\pos acc -> if length pos == 1 then head pos : acc else acc) [] (positionMap k)

prop_suboptimalRuleIntervalKeepBalancedGenomes :: Property
prop_suboptimalRuleIntervalKeepBalancedGenomes = suboptimalRuleKeepBalancedGenomes suboptimalRuleInterval

prop_suboptimalRulePairsKeepBalancedGenomes :: Property
prop_suboptimalRulePairsKeepBalancedGenomes = suboptimalRuleKeepBalancedGenomes suboptimalRulePairs

suboptimalRuleKeepBalancedGenomes :: (RigidFlexibleReverseMatcher GenesIRsR GenesIRsF -> GenesIRsR -> GenesIRsF -> (GenesIRsR, GenesIRsF)) -> Property
suboptimalRuleKeepBalancedGenomes rule =
  property $ do
    g <- forAll (genRGenome 100)
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h <- forAll $ rearrangeAndFlexibilizeGenome k g
    (g', h') <- forAll . return $ rule RFRM g h
    assert $ areBalanced RFRM g' h'

prop_bpsToBlockDelsIsomorphism :: Property
prop_bpsToBlockDelsIsomorphism = property $ do
  (GW g) <- forAll (genGenome 100)
  k <- forAll $ Gen.int (Range.linear 0 (size g - 2))
  bps_ <- forAll . fmap (take k) . Gen.shuffle $ [1 .. mkIdx (size g - 1)]
  let bps = EnumSet.fromList bps_
  bps === (blockDelsToBps . bpsToBlockDels g $ bps)

prop_equalGenomesHaveOnlyPerfectComponetsOnBMG :: Property
prop_equalGenomesHaveOnlyPerfectComponetsOnBMG =
  property $ do
    g <- forAll (genRGenome 100)
    let pg = trivialPartition g
    let (comps, _) = getConnectedComponents (mkBlockMatchGraph RRRM pg pg)
    assert $ countPerfect comps == IntMap.size comps

prop_greedyPartitionProduceValidCorrespondence :: Property
prop_greedyPartitionProduceValidCorrespondence =
  partitionProduceValidCorrespondence (greedyPartition False) 100

prop_greedyPartitionSingProduceValidCorrespondence :: Property
prop_greedyPartitionSingProduceValidCorrespondence =
  partitionProduceValidCorrespondence (greedyPartition True) 100

prop_soarPartitionProduceValidCorrespondence :: Property
prop_soarPartitionProduceValidCorrespondence =
  partitionProduceValidCorrespondence (soarPartition False) 100

prop_approxPartitionProduceValidCorrespondence :: Property
prop_approxPartitionProduceValidCorrespondence =
  partitionProduceValidCorrespondence approxPartition 50

partitionProduceValidCorrespondence :: (RigidFlexibleReverseMatcher GenesIRsR GenesIRsF -> GenesIRsR -> GenesIRsF -> CommonPartition GenesIRsR GenesIRsF) -> Int -> Property
partitionProduceValidCorrespondence partAlg size_lim =
  property $ do
    g <- forAll (genRGenome size_lim)
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h <- forAll $ rearrangeAndFlexibilizeGenome k g
    part <- forAll . return $ partAlg RFRM g h
    (CGP pg ph) <- forAll . return $ combine RFRM part
    assert $ checkCommon RFRM pg ph

prop_fptPartitionProduceValidCorrespondence :: Property
prop_fptPartitionProduceValidCorrespondence =
  property $ do
    g <- forAll (genRGenome 20)
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h <- forAll $ rearrangeAndFlexibilizeGenome k g
    (part,_) <- fmap (first fromJust) . evalIO $ fptPartition 100000000 RFRM g h
    (CGP pg ph) <- forAll . return $ combine RFRM part
    assert $ checkCommon RFRM pg ph

tests :: IO Bool
tests = checkSequential $$(discover)
-- tests = check prop_fptPartitionProduceValidCorrespondence
-- tests = do
  -- recheck (Size 0) (Seed 12439291806607096571 15074968846026035773) prop_equalGenomesHaveOnlyPerfectComponetsOnBMG
  -- return True
