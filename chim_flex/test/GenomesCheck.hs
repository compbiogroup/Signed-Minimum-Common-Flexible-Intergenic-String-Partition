{-# LANGUAGE TemplateHaskell #-}

module GenomesCheck (tests, genRGenome, rearrangeGenome) where

import Control.Monad.Random (MonadRandom, getRandomR, replicateM)
import Genomes
import Hedgehog
import qualified Hedgehog.Gen as Gen
import qualified Hedgehog.Range as Range
import LocalBase

genGene :: Gen Int
genGene = Gen.int (Range.linear 1 100)

genIR :: Gen Int
genIR = Gen.int (Range.linear 0 100)

-- | Generate an empty genome
genRGenome :: Gen GenesIRsR
genRGenome = genRGenomeWithSign True

genRGenomeWithSign :: Bool -> Gen GenesIRsR
genRGenomeWithSign signed = do
  n <- Gen.int (Range.linear 1 100)
  coins <- Gen.list (Range.singleton n) Gen.bool
  genes <- (if signed
          then zipWith swaps coins
          else id)
      <$> Gen.list (Range.singleton n) genGene
  irs <- Gen.list (Range.singleton $ n - 1) genIR
  return $ mkRGenome genes irs
  where
    swaps b v = if b then v else -v
    
rearrangeGenome :: (MonadRandom mon, RigidIntergenicGenome g) => Int -> g -> mon g
rearrangeGenome k g =
  do
    revs <- replicateM k rev
    return $ foldr (\(i,j) g' -> intergenicReversal (mkIdx i) (mkIdx j) 0 0 g') g revs
  where
    rev = do
      i <- getRandomR (2::Int, size g - 1)
      j <- getRandomR (i, size g - 1)
      return (i,j)

tests :: IO Bool
tests = checkSequential $$(discover)
