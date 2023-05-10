{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE ImportQualifiedPost #-}

module PartitionCheck (tests) where

import Control.Monad.Random (evalRandIO)
import Data.IntSet qualified as IntSet
import Data.List qualified as List
import Data.Maybe (fromJust)
import Genomes
import GenomesCheck (GenomeWrapper (..), genGenome, genRGenome, rearrangeGenome)
import Hedgehog
import Hedgehog.Gen qualified as Gen
import Hedgehog.Range qualified as Range
import LocalBase
import Partition

prop_commonPrefixFromBothAreEqual :: Property
prop_commonPrefixFromBothAreEqual =
  property $ do
    g <- forAll genRGenome
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h <- evalIO . evalRandIO $ rearrangeGenome k g
    idx_g <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
    idx_h <- forAll . fmap head . Gen.shuffle . fromJust $ geneMapLookup (getGene idx_g g) (positionMap h)
    let ((g_beg, g_end), (h_beg, h_end)) = fromJust $ commonPrefix Nothing RRRM IntSet.empty IntSet.empty g h idx_g idx_h
    subG <- forAll . return $ subGenome g_beg g_end g
    subH <- forAll . return $ subGenome h_beg h_end h
    assert $ isMatch RRRM subG subH

prop_longestSubstringFromBothAreEqual :: Property
prop_longestSubstringFromBothAreEqual =
  property $ do
    g <- forAll genRGenome
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h_ <- evalIO . evalRandIO $ rearrangeGenome k g
    h <- forAll . return $ h_
    let (((g_beg, g_end), (h_beg, h_end)), _, _) = fromJust $ longestSubstring Nothing RRRM IntSet.empty IntSet.empty g h
    let subG = subGenome g_beg g_end g
        subH = subGenome h_beg h_end h
    assert $ isMatch RRRM subG subH

prop_suboptimalRuleIntervalTurnsReplicasIntoSingletons :: Property
prop_suboptimalRuleIntervalTurnsReplicasIntoSingletons = suboptimalRuleTurnsReplicasIntoSingletons suboptimalRuleInterval

prop_suboptimalRulePairsTurnsReplicasIntoSingletons :: Property
prop_suboptimalRulePairsTurnsReplicasIntoSingletons = suboptimalRuleTurnsReplicasIntoSingletons suboptimalRulePairs

suboptimalRuleTurnsReplicasIntoSingletons :: (RigidRigidReverseMatcher GenesIRsR GenesIRsR -> GenesIRsR -> GenesIRsR -> (GenesIRsR, GenesIRsR)) -> Property
suboptimalRuleTurnsReplicasIntoSingletons rule =
  property $ do
    g <- forAll genRGenome
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h_ <- evalIO . evalRandIO $ rearrangeGenome k g
    h <- forAll . return $ h_
    (g', h') <- forAll . return $ rule RRRM g h
    let sigG = getSingletons g
        sigH = getSingletons h
        sigG' = getSingletons g'
        sigH' = getSingletons h'
    assert $ sigG `List.isSubsequenceOf` sigG'
    assert $ sigH `List.isSubsequenceOf` sigH'
  where
    getSingletons k = List.sort $ foldr (\pos acc -> if length pos == 1 then head pos : acc else acc) [] (positionMap k)

prop_bpsToBlocksIsomorphism :: Property
prop_bpsToBlocksIsomorphism = property $ do
  (GW g) <- forAll genGenome
  k <- forAll $ Gen.int (Range.linear 0 (size g - 2))
  bps_ <- forAll . fmap (take k) . Gen.shuffle $ [1 .. mkIdx (size g - 1)]
  let bps = (0) : bps_ ++ [mkIdx (size g)]
  bps === (blocksToBps . bpsToBlocks g $ bps)

tests :: IO Bool
tests = checkSequential $$(discover)
