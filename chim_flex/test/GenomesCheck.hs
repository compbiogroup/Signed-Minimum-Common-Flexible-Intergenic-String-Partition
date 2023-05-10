{-# LANGUAGE ExistentialQuantification #-}
{-# LANGUAGE TemplateHaskell #-}

module GenomesCheck (tests, genGenome, genRGenome, rearrangeGenome, GenomeWrapper (..)) where

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

data GenomeWrapper = forall g. (Genome g, Show g) => GW g

instance Show GenomeWrapper where
  show (GW g) = show g

-- Genomes will have size at least 3 (extensions plus at least one gene)
genGenome :: Gen GenomeWrapper
genGenome = do
  flex <- Gen.bool
  if flex
    then GW <$> genFGenome
    else GW <$> genRGenome

-- | Generate an empty rigid genome
genRGenome :: Gen GenesIRsR
genRGenome = genRGenomeWithSign Signed

-- | Generate an empty flexible genome
genFGenome :: Gen GenesIRsF
genFGenome = genFGenomeWithSign Signed

genRGenomeWithSign :: Sign -> Gen GenesIRsR
genRGenomeWithSign sign = do
  n <- Gen.int (Range.linear 3 100)
  coins <- Gen.list (Range.singleton n) Gen.bool
  genes <-
    ( case sign of
        Signed -> zipWith swaps coins
        Unsigned -> id
      )
      <$> Gen.list (Range.singleton $ n - 2) genGene
  irs <- Gen.list (Range.singleton $ n - 1) genIR
  return $ mkRGenome True sign genes irs
  where
    swaps b v = if b then v else -v

genFGenomeWithSign :: Sign -> Gen GenesIRsF
genFGenomeWithSign sign = do
  n <- Gen.int (Range.linear 3 100)
  coins <- Gen.list (Range.singleton n) Gen.bool
  genes <-
    ( case sign of
        Signed -> zipWith swaps coins
        Unsigned -> id
      )
      <$> Gen.list (Range.singleton $ n - 2) genGene
  irs_low <- Gen.list (Range.singleton $ n - 1) genIR
  porc <- Gen.list (Range.singleton $ n - 1) $ Gen.int (Range.linear 1 50)
  let irs_high = zipWith (\ir p -> ir + p * ir `div` 100) irs_low porc
  let irs = zip irs_low irs_high
  return $ mkFGenome True sign genes irs
  where
    swaps b v = if b then v else -v

rearrangeGenome :: (MonadRandom mon, RigidIntergenicGenome g) => Int -> g -> mon g
rearrangeGenome k g =
  do
    revs <- replicateM k rev
    return $ foldr (\(i, j) g' -> intergenicReversal (mkIdx i) (mkIdx j) 0 0 g') g revs
  where
    rev = do
      i <- getRandomR (2 :: Int, size g - 1)
      j <- getRandomR (i, size g - 1)
      return (i, j)

prop_getGeneIsGene :: Property
prop_getGeneIsGene = property $ do
  (GW g) <- forAll genGenome
  idx <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  assert $ getGene idx g `isGene` g

prop_invGeneIsGene :: Property
prop_invGeneIsGene = property $ do
  (GW g) <- forAll genGenome
  idx <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  assert $ invGene g (getGene idx g) `isGene` g

prop_invGeneIsAutoinverse :: Property
prop_invGeneIsAutoinverse = property $ do
  (GW g) <- forAll genGenome
  idx <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  let gene = getGene idx g
  invGene g (invGene g gene) === gene

prop_setGenesAreGenesInIndices :: Property
prop_setGenesAreGenesInIndices = property $ do
  (GW g) <- forAll genGenome
  genes <- forAll $ Gen.list (Range.linear 1 (size g - 2)) (mkGene <$> Gen.int (Range.linear 1 200))
  indices <- forAll $ take (length genes) <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
  let g' = setGenes indices genes g
  map (`getGene` g') indices === genes

prop_makeSingletonsAddSingletonsInIndices :: Property
prop_makeSingletonsAddSingletonsInIndices = property $ do
  (GW g) <- forAll genGenome
  k <- forAll $ Gen.int (Range.linear 1 (min 5 (size g - 2)))
  indices <- forAll $ take k <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
  let (g', genes) = makeSingletons g indices g
  let posMap = positionMap g'
      testSingleton gene =
        case geneMapLookup gene posMap of
          Nothing -> False
          Just pos -> length pos == 1
  map (`getGene` g') indices === genes
  assert $ all testSingleton genes

prop_alphabetAreGenes :: Property
prop_alphabetAreGenes = property $ do
  (GW g) <- forAll genGenome
  let alp = alphabet g
  assert $ all (`isGene` g) alp

prop_newGeneIsNotGeneAftersetGenes :: Property
prop_newGeneIsNotGeneAftersetGenes =
  property $ do
    (GW g) <- forAll genGenome
    genes <- forAll $ Gen.list (Range.linear 1 (size g - 2)) (mkGene <$> Gen.int (Range.linear 1 200))
    indices <- forAll $ take (length genes) <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
    let g' = setGenes indices genes g
    assert $ not (getNewGene g' `isGene` g')

prop_newGeneIsNotGeneAfterMakeSingletons :: Property
prop_newGeneIsNotGeneAfterMakeSingletons = property $ do
  (GW g) <- forAll genGenome
  k <- forAll $ Gen.int (Range.linear 1 (min 5 (size g - 2)))
  indices <- forAll $ take k <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
  let (g', genes) = makeSingletons g indices g
  assert $ not (getNewGene g' `isGene` g')

tests :: IO Bool
tests = checkSequential $$(discover)