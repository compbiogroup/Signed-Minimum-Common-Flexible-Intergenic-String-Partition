{-# LANGUAGE ExistentialQuantification #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE TemplateHaskell #-}

module GenomesCheck (tests, genGenome, genRGenome, rearrangeGenome, rearrangeAndFlexibilizeGenome, GenomeWrapper (..)) where

import Control.Monad.Random (replicateM)
import Genomes
import Hedgehog
import Hedgehog.Gen qualified as Gen
import Hedgehog.Range qualified as Range
import LocalBase

genGene :: Gen Int
genGene = Gen.int (Range.linear 1 100)

genIR :: Gen Int
genIR = Gen.int (Range.linear 0 100)

data GenomeWrapper = forall g. (Genome g, Show g) => GW g

instance Show GenomeWrapper where
  show (GW g) = show g

-- Genomes will have size at least 3 (extensions plus at least one gene)
genGenome :: Int -> Gen GenomeWrapper
genGenome size_lim = do
  flex <- Gen.bool
  if flex
    then GW <$> genFGenome size_lim
    else GW <$> genRGenome size_lim

-- | Generate an empty rigid genome
genRGenome :: Int -> Gen GenesIRsR
genRGenome = genRGenomeWithSign Signed

-- | Generate an empty flexible genome
genFGenome :: Int -> Gen GenesIRsF
genFGenome = genFGenomeWithSign Signed

genRGenomeWithSign :: Sign -> Int -> Gen GenesIRsR
genRGenomeWithSign sign size_lim = do
  n <- Gen.int (Range.linear 3 size_lim)
  coins <- Gen.list (Range.singleton n) Gen.bool
  genes <-
    ( case sign of
        Signed -> zipWith swaps coins
        Unsigned -> id
      )
      <$> Gen.list (Range.singleton $ n - 2) genGene
  irs <- Gen.list (Range.singleton $ n - 1) genIR
  return $ mkRGenome Linear True sign genes irs
  where
    swaps b v = if b then v else -v

genFGenomeWithSign :: Sign -> Int -> Gen GenesIRsF
genFGenomeWithSign sign size_lim = do
  n <- Gen.int (Range.linear 3 size_lim)
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
  return $ mkFGenome Linear True sign genes irs
  where
    swaps b v = if b then v else -v

rearrangeGenome :: (RigidIntergenicChromosome g) => Int -> g -> Gen g
rearrangeGenome = applyReversals

rearrangeAndFlexibilizeGenome :: Int -> GenesIRsR -> Gen GenesIRsF
rearrangeAndFlexibilizeGenome k g = do
  g' <- rearrangeGenome k g
  f <- Gen.int (Range.linear 0 100)
  return (flexibilize f g')

applyReversals :: (RigidIntergenicChromosome g) => Int -> g -> Gen g
applyReversals k g =
  if size g <= 4
    then return g
    else do
      revs <- replicateM k rev
      return $ foldr (\(i, j) g' -> intergenicReversal (mkIdx i) (mkIdx j) 0 0 g') g revs
  where
    rev = do
      i <- Gen.int (Range.linear (2 :: Int) (size g - 2))
      j <- Gen.int (Range.linear (i + 1) (size g - 1))
      return (i, j)

prop_getGeneIsGene :: Property
prop_getGeneIsGene = property $ do
  (GW g) <- forAll (genGenome 100)
  idx <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  assert $ getGene idx g `isGene` g

prop_invGeneIsGene :: Property
prop_invGeneIsGene = property $ do
  (GW g) <- forAll (genGenome 100)
  idx <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  assert $ invGene g (getGene idx g) `isGene` g

prop_invGeneIsAutoinverse :: Property
prop_invGeneIsAutoinverse = property $ do
  (GW g) <- forAll (genGenome 100)
  idx <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  let gene = getGene idx g
  invGene g (invGene g gene) === gene

prop_setGenesAreGenesInIndices :: Property
prop_setGenesAreGenesInIndices = property $ do
  (GW g) <- forAll (genGenome 100)
  genes <- forAll $ Gen.list (Range.linear 1 (size g - 2)) (mkGene <$> Gen.int (Range.linear 1 200))
  indices <- forAll $ take (length genes) <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
  let g' = setGenes indices genes g
  map (`getGene` g') indices === genes

prop_makeSingletonsAddSingletonsInIndices :: Property
prop_makeSingletonsAddSingletonsInIndices = property $ do
  (GW g) <- forAll (genGenome 100)
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
  (GW g) <- forAll (genGenome 100)
  let alp = alphabet g
  assert $ all (`isGene` g) alp

prop_newGeneIsNotGeneAftersetGenes :: Property
prop_newGeneIsNotGeneAftersetGenes =
  property $ do
    (GW g) <- forAll (genGenome 100)
    genes <- forAll $ Gen.list (Range.linear 1 (size g - 2)) (mkGene <$> Gen.int (Range.linear 1 200))
    indices <- forAll $ take (length genes) <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
    let g' = setGenes indices genes g
    assert $ not (getNewGene g' `isGene` g')

prop_newGeneIsNotGeneAfterMakeSingletons :: Property
prop_newGeneIsNotGeneAfterMakeSingletons = property $ do
  (GW g) <- forAll (genGenome 100)
  k <- forAll $ Gen.int (Range.linear 1 (min 5 (size g - 2)))
  indices <- forAll $ take k <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
  let (g', _) = makeSingletons g indices g
  assert $ not (getNewGene g' `isGene` g')

prop_reversalsKeepBalanced :: Property
prop_reversalsKeepBalanced = property $ do
  g <- forAll (genRGenome 100)
  k <- forAll $ Gen.int (Range.linear 0 (size g))
  h <- forAll $ applyReversals k g
  assert $ areBalanced RRRM g h

prop_occurrenceWithGeneMapImplementationWorks :: Property
prop_occurrenceWithGeneMapImplementationWorks = property $ do
  g <- forAll (genRGenome 100)
  i <- forAll (Gen.int (Range.linear 1 (size g)))
  let a = canonicOri (getGene (mkIdx i) g)
  occurrence (positionMap g) a === (sum . fmap (\i' -> if abs (getGene (mkIdx i') g) == a then 1 else 0) $ [1 .. size g])

prop_joinIsInverseOfCut :: Property
prop_joinIsInverseOfCut = property $ do
  g <- toMC <$> forAll (genRGenome 100)
  chr_i <- mkIdx <$> forAll (Gen.int (Range.linear (1 :: Int) (numChromosomes g)))
  i <- mkIdx <$> forAll (Gen.int (Range.linear (1 :: Int) (size $ getChromosome chr_i g)))
  (join chr_i (chr_i + 1) False False False . cut (chr_i,i) 0 $ g) === g

tests :: IO Bool
tests = checkSequential $$(discover)
