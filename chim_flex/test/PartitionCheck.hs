{-# LANGUAGE TemplateHaskell #-}

module PartitionCheck (tests) where

import Control.Monad.Random (evalRandIO)
import qualified Data.IntSet as IntSet
import Data.Maybe (fromJust)
import Genomes (Genome (size, subGenome), Matcher (isMatch), RigidRigidReverseMatcher (RRRM))
import GenomesCheck (genRGenome, rearrangeGenome)
import Hedgehog
  ( Property,
    assert,
    checkSequential,
    discover,
    evalIO,
    forAll,
    property,
  )
import qualified Hedgehog.Gen as Gen
import qualified Hedgehog.Range as Range
import LocalBase
import Partition

prop_longestSubstringFromBothAreEqual :: Property
prop_longestSubstringFromBothAreEqual =
  property $ do
    g <- forAll genRGenome
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h <- evalIO . evalRandIO $ rearrangeGenome k g
    let (((g_beg, g_end), (h_beg, h_end)), _, _) = fromJust $ longestSubstring Nothing RRRM IntSet.empty IntSet.empty g h
    assert $ isMatch RRRM (subGenome g_beg g_end g) (subGenome h_beg h_end h)

tests :: IO Bool
tests = checkSequential $$(discover)
