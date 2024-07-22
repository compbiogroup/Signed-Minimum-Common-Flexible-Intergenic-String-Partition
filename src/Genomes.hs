{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE ExplicitForAll #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

-- |
-- Module      : Genomes
-- Description : Representation of a Genome. The representation may include information
-- regarding intergenic regions and gene orientation.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module Genomes
  ( Genome (..),
    Chromosome (..),
    MultiChromosome (..),
    IntergenicChromosome (..),
    RigidIntergenicChromosome (..),
    RigidIntergenicMultiChromGenome (..),
    Gene,
    mkGene,
    geneToInt,
    IR,
    irToInt,
    Idx,
    Sign (..),
    ChromType(..),
    mkIdx,
    idxDist,
    idxToInt,
    incIdx,
    decIdx,
    occurrence,
    GeneMap,
    positionMap,
    geneMapLookup,
    geneMapAdjust,
    occurrenceMax,
    singletonOnBoth,
    MultiGIs,
    GenesIRsR,
    GenesIRsF,
    mkMCRGenome,
    mkRGenome,
    mkFGenome,
    mkRGenome0,
    mkFGenome0,
    writeIR,
    readRGenome,
    writeRGenome,
    readFGenome,
    writeFGenome,
    writeMultiC,
    flexibilize,
    flexibilizeMC,
    Matcher (..),
    FlipMatcher (..),
    RigidRigidDirectMatcher (..),
    RigidFlexibleDirectMatcher (..),
    RigidRigidReverseMatcher (..),
    RigidFlexibleReverseMatcher (..),
    isCompatibleWithSubGenome,
    randomGenome,
    randomGenomeWithReplicas,
    shuffleGenome,
  )
where

import Control.Exception (assert)
import Control.Monad.Random (MonadRandom, getRandomRs, getRandoms)
import Control.Monad.ST (ST)
import Data.ByteString.Builder (intDec, toLazyByteString)
import Data.ByteString.Char8 qualified as BS
import Data.ByteString.Lazy qualified as LBS
import Data.Coerce (coerce)
import Data.Hashable (Hashable)
import Data.IntMap (IntMap)
import Data.IntMap qualified as IntMap
import Data.List (find, intersperse)
import Data.List qualified as List
import Data.Maybe (fromJust)
import Data.Set (Set)
import Data.Set qualified as Set
import Data.Vector (Vector, (!))
import Data.Vector qualified as Vec
import Data.Vector.Mutable qualified as MVec
import LocalBase
import System.Random (Random)
import System.Random.Shuffle (shuffleM)

newtype Gene = Gene Int deriving newtype (Eq, Show, Read, Hashable, Ord, Num, Bounded, Enum, Random)

mkGene :: Int -> Gene
mkGene = Gene

geneToInt :: Gene -> Int
geneToInt (Gene i) = i

data IR = R Int | F Int Int deriving (Eq, Ord)

irToInt :: IR -> Int
irToInt ir = case ir of
  R i -> i
  F _ _ -> error patternError

sumIR :: IR -> IR -> IR
sumIR (R i) (R j) = R (i+j)
sumIR _ _ = error patternError

newtype Idx = Idx Int deriving newtype (Eq, Show, Read, Hashable, Ord, Num, Bounded, Enum, Integral, Real, Random)

mkIdx :: Int -> Idx
mkIdx = Idx

idxToInt :: Idx -> Int
idxToInt (Idx idx) = idx

idxDist :: Idx -> Idx -> Int
-- Size of the interval containing the two indices, including the indices
idxDist (Idx idx1) (Idx idx2) = idx2 - idx1 + 1

incIdx :: Idx -> Idx
incIdx (Idx idx) = Idx (idx + 1)

decIdx :: Idx -> Idx
decIdx (Idx idx) = Idx (idx - 1)

instance Show IR where
  show (R x) = show x
  show (F x y) = show x ++ ":" ++ show y

instance Orientable Gene where
  getOri a = if a >= 0 then LR else RL
  invOri a = -a
  hasOri _ = True

data Sign = Signed | Unsigned deriving (Eq, Show, Enum)

data ChromType = Circular | Linear deriving (Eq, Show, Enum)

-- Representation with a list of genes and a list of intergenic regions
-- (only ints for rigid, or intervals for flexible)
-- We assume that linear chromosomes start and end with genes (which are
-- different from all other genes unless you use turnLinear) and the
-- intergenic region between the first and last gene in the list of a
-- circular chromosomes is the last one in the intergenic region list.
data GenesIRs = GenesIRs ChromType Sign (Vector Gene) (Vector IR) Gene

instance Show GenesIRs where
  show (GenesIRs Linear _ genes irs _) =
    unwords . (("(" ++ head str_s ++ ")") :) $
      zipWith (\ir a -> "- " ++ ir ++ " - (" ++ a ++ ")") str_i (tail str_s)
    where
      str_s = Vec.toList $ (\i -> if i == maxBound then "inf" else show i) <$> genes
      str_i = Vec.toList $ show <$> irs
  show (GenesIRs Circular _ genes irs _) =
    unwords $ zipWith (\ir a -> "(" ++ a ++ ") - " ++ ir ++ " - ") str_i str_s
    where
      str_s = Vec.toList $ (\i -> if i == maxBound then "inf" else show i) <$> genes
      str_i = Vec.toList $ show <$> irs

-- Representation with multiple lists of genes and multiple lists of intergenic regions
-- (only ints for rigid, or intervals for flexible). Each list represent a linear or circular chromosome.
-- make sure the genes on the extremities of linear chromosomes are always in the caps list.
data MultiGIs c = MultiGIs {caps :: [Gene], chromosomes :: [c]}

instance (Show c) => Show (MultiGIs c) where
  show (MultiGIs _ chroms) = unwords . intersperse " | " . map show $ chroms

-- TODO: make Idx a TypeFamily of Genome
class Genome g where
  isGene :: Gene -> g -> Bool

  -- Number of genes
  size :: g -> Int

  -- get gene at position i (index starts in 1), requires 1 <= i <= size g
  getGene :: Idx -> g -> Gene

  -- invert gene according to the characteristics of the genome
  -- (takes into account if the genome has known orientation)
  invGene :: g -> Gene -> Gene

  -- set genes at indices requires 1 < i < size g,
  -- for each index i.
  -- The signs must be specify on the new genes,
  -- because the original sings will be overwritten
  setGenes :: [Idx] -> [Gene] -> g -> g

  -- for indices i j, it gets the subgenome
  -- starting at gene position i
  -- and ending at gene position j
  -- requires 1 <= i < j <= size g
  subGenome :: Idx -> Idx -> g -> g

  -- Return set with all the gene values
  alphabet :: g -> Set Gene

  -- Get a gene that does not appear on the genome (using Int as the
  -- underline representation, this value should be bigger than all
  -- other gene values)
  getNewGene :: g -> Gene

  -- Replace genes in indices with a new singleton gene
  -- (does not have any replica in the genome),
  -- Returns the genomes and the singleton that replace each index.
  -- Also ensures that they are singletons in g2.
  -- Requires 1 < i < size g, for each index i.
  makeSingletons :: (Genome g2) => g2 -> [Idx] -> g -> (g, [Gene])

class (Genome g) => Chromosome g where
  isLinear :: g -> Bool
  isCircular :: g -> Bool
  isCircular = not . isLinear
  turnCircular :: g -> g

  -- Turn into a linear chromossome by cutting after the position indicated by the index
  turnLinear :: Idx -> g -> g

  -- Combine two linear chromosomes into one putting a intergenic region between them.
  combine :: g -> IR -> g -> g

class (Genome g) => MultiChromosome g where
  type Chrom g
  getChromosome :: Idx -> g -> Chrom g
  numChromosomes :: g -> Int

class (Chromosome g) => IntergenicChromosome g where
  intergenicFullReversal :: g -> g
  getIR :: Idx -> g -> IR

class (IntergenicChromosome g) => RigidIntergenicChromosome g where
  intergenicReversal :: Idx -> Idx -> Int -> Int -> g -> g
  intergenicTransposition :: Idx -> Idx -> Idx -> Int -> Int -> Int -> g -> g
  intergenicInsertion :: Idx -> g -> g -> g
  intergenicDeletion :: Idx -> Idx -> Int -> g -> g

class (MultiChromosome g) => RigidIntergenicMultiChromGenome g where
  -- Cut between (idxi,idxi+1) = (a,b) and (idxj,idxj+1) = (c,d) followed by join of (a,c) and (b,d) or (a,d) and (b,c)
  intergenicDCJ :: (Idx, Idx) -> (Idx, Idx) -> Bool -> Int -> Int -> g -> g

  -- Cut the chromossome indicated by the first index after the position indicated
  -- by the second index, the cut is done after the given number of nucleotides.
  -- Add caps accordingly. Circular chromosomes are turned linear by this operation.
  cut :: (Idx, Idx) -> Int -> g -> g

  -- Join the chromosomes in the given indices (remove caps accordingly).
  -- The first bool indicates the direction (True to move the first and False to move the second).
  -- The second and third bools indicate whether the correspondent chromosome should be inverted.
  join :: Idx -> Idx -> Bool -> Bool -> Bool -> g -> g

  -- Turn chromosome at position into a circular chromosome
  turnCircularAtIndex :: Idx -> g -> g

occurrence :: GeneMap [Idx] -> Gene -> Int
occurrence geneMap a = maybe 0 length (geneMapLookup a geneMap)

newtype GeneMap a = GM (IntMap a) deriving newtype (Functor, Foldable)

positionMap :: (Genome g) => g -> GeneMap [Idx]
positionMap g = GM . IntMap.fromListWith (++) . map (\i -> (geneToInt . abs $ getGene i g, [i])) $ [1 .. mkIdx (size g)]

geneMapLookup :: Gene -> GeneMap a -> Maybe a
geneMapLookup a (GM gm) = IntMap.lookup (geneToInt $ abs a) gm

geneMapAdjust :: (a -> a) -> Gene -> GeneMap a -> GeneMap a
geneMapAdjust f a (GM gm) = GM $ IntMap.adjust f (geneToInt $ abs a) gm

occurrenceMax :: (Genome g) => g -> Int
occurrenceMax = maximum . fmap length . positionMap

singletonOnBoth :: GeneMap [Idx] -> GeneMap [Idx] -> Gene -> Bool
-- ^ Check if the gene is a singleton in both genomes
singletonOnBoth posMapG posMapH gene =
  (case geneMapLookup gene posMapG of Nothing -> False; Just pos -> length pos == 1)
    && (case geneMapLookup gene posMapH of Nothing -> False; Just pos -> length pos == 1)

instance Genome GenesIRs where
  isGene gene (GenesIRs _ sign genes _ _) =
    gene `elem` genes
      || case sign of
        Signed -> (-gene) `elem` genes
        Unsigned -> False
  size (GenesIRs _ _ genes _ _) = Vec.length genes
  getGene idx (GenesIRs _ _ genes _ _) = genes Vec.! (idxToInt idx - 1)

  invGene (GenesIRs _ sign _ _ _) a =
    case sign of Signed -> -a; Unsigned -> a

  setGenes idxs new_genes (GenesIRs chtype sign genes irs new) = GenesIRs chtype sign genes' irs new'
    where
      genes' = Vec.modify update genes
      update :: MVec.MVector s Gene -> ST s ()
      update v =
        mapM_ (\(idx, gene) -> MVec.write v (idxToInt idx - 1) gene) $
          zip idxs new_genes
      new' = foldr (\gene n -> if abs gene >= n then abs gene + 1 else n) new new_genes

  subGenome idxi idxj (GenesIRs _ sign genes irs new) = GenesIRs Linear sign (Vec.slice i n genes) (Vec.slice i (n - 1) irs) new
    where
      i = idxToInt idxi - 1
      j = idxToInt idxj
      n = j - i
  alphabet (GenesIRs _ _ genes _ _) = Set.fromList . Vec.toList . fmap abs $ genes

  -- In this implementation the new gene is bigger than all other genes except for maxBound
  getNewGene (GenesIRs _ _ _ _ new) = new

  makeSingletons h idxs g@(GenesIRs _ _ _ _ new_g) = (g', new_genes)
    where
      g' = setGenes idxs new_genes g
      new_h = getNewGene h
      new_genes = take (length idxs) [max new_g new_h ..]

instance Chromosome GenesIRs where
  -- turn into a circular chromosome. The first and last genes (caps) are deleted and the first
  -- and last IRs are combined.
  turnCircular (GenesIRs Linear sign genes irs new) = GenesIRs Circular sign (Vec.init . Vec.tail $ genes) (merge_ir irs) new
    where
      merge_ir v =
        let n = Vec.length v in
          Vec.tail $ v Vec.// [(n - 1, sumIR (v Vec.! 0) (v Vec.! (n - 1)))]
  turnCircular g = g
  turnLinear idx (GenesIRs Circular sign genes irs new) = GenesIRs Linear sign genes' irs' new
    where
      genes' = Vec.fromList . rotateL (idxToInt idx) . Vec.toList $ genes
      irs' = Vec.fromList . tail . rotateL (idxToInt idx - 1) . Vec.toList $ irs
  turnLinear _ (GenesIRs Linear _ _ _ _) = error patternError
  combine (GenesIRs Circular _ _ _ _) _ _ = error logicError
  combine _ _ (GenesIRs Circular _ _ _ _) = error logicError
  combine (GenesIRs Linear sign1 genes1 irs1 new1) ir (GenesIRs Linear sign2 genes2 irs2 new2) =
    if sign1 /= sign2 then error logicError else GenesIRs Linear sign1 genes' irs' new'
    where
      genes' = Vec.concat [genes1, genes2]
      irs' = Vec.concat [irs1, Vec.singleton ir, irs2]
      new' = max new1 new2
  isLinear (GenesIRs Linear _ _ _ _) = True
  isLinear (GenesIRs Circular _ _ _ _) = False

instance Genome c => Genome (MultiGIs c) where
  isGene gene = any (isGene gene) . chromosomes
  size = sum . map size . chromosomes
  getGene i = findGene i . chromosomes
    where
      findGene idx (g : gs) =
        let s = size g
         in if idx <= mkIdx s then getGene idx g else findGene (idx - mkIdx s) gs
      findGene _ [] = error indexError
  invGene g = invGene (head . chromosomes $ g)
  setGenes idxs genes mc = foldl setGene mc (zip idxs genes)
    where
      setGene genome (i, gene) = MultiGIs (caps mc) . findChrom [] i . chromosomes $ genome
        where
          findChrom gs_old idx (g : gs) =
            let s = size g
             in if idx <= mkIdx s then reverse gs_old ++ (setGenes [idx] [gene] g : gs) else findChrom (g : gs_old) (idx - mkIdx s) gs
          findChrom _ _ [] = error indexError
  subGenome idxi0 idxj0 mc = MultiGIs (caps mc) . findChrom [] idxi0 idxj0 . chromosomes $ mc
    where
      findChrom gs_old idxi idxj (g : gs) =
        let s = size g
         in if
                | idxj <= mkIdx s -> reverse (subGenome idxi idxj g : gs_old)
                | idxi == 0 -> findChrom (g : gs_old) idxi (idxj - mkIdx s) gs
                | idxi <= mkIdx s -> findChrom [subGenome idxi (mkIdx s) g] 0 (idxj - mkIdx s) gs
                | otherwise -> findChrom [] (idxi - mkIdx s) (idxj - mkIdx s) gs
      findChrom _ _ _ [] = error indexError
  alphabet = Set.unions . map alphabet . chromosomes
  getNewGene = maximum . map getNewGene . chromosomes
  makeSingletons h idxs g = (g', new_genes)
    where
      g' = setGenes idxs new_genes g
      new_g = getNewGene g
      new_h = getNewGene h
      new_genes = take (length idxs) [max new_g new_h ..]

instance IntergenicChromosome GenesIRs where
  getIR idx (GenesIRs _ _ _ irs _) = irs Vec.! (idxToInt idx - 1)
  intergenicFullReversal (GenesIRs Linear sign vs vi new) = GenesIRs Linear sign vs' vi' new
    where
      vs' = case sign of { Unsigned -> id; Signed -> fmap negate } $ Vec.reverse vs
      vi' = Vec.reverse vi
  intergenicFullReversal (GenesIRs Circular sign vs vi new) = GenesIRs Circular sign vs' vi new
    where
      vs' = case sign of { Unsigned -> id; Signed -> fmap negate } $ vs

instance Orientable GenesIRs where
  getOri g@(GenesIRs _ _ genes irs _) = if (genes, irs) <= (genes', irs') then LR else RL
    where
      (GenesIRs _ _ genes' irs' _) = invOri g
  invOri = intergenicFullReversal
  hasOri (GenesIRs _ sign _ _ _) = sign == Signed

newtype GenesIRsR = GLR GenesIRs deriving newtype (Show, Genome, Chromosome, IntergenicChromosome, Orientable)

newtype GenesIRsF = GLF GenesIRs deriving newtype (Show, Genome, Chromosome, IntergenicChromosome, Orientable)

instance Eq GenesIRsR where
  (GLR (GenesIRs Linear sign1 genes1 irs1 _)) == (GLR (GenesIRs Linear sign2 genes2 irs2 _)) =
    sign1 == sign2 && genes1 == genes2 && irs1 == irs2
  (GLR (GenesIRs Circular _ _ _ _)) == (GLR (GenesIRs Circular _ _ _ _)) = undefined
  (GLR (GenesIRs Linear _ _ _ _)) == (GLR (GenesIRs Circular _ _ _ _)) = False
  (GLR (GenesIRs Circular _ _ _ _)) == (GLR (GenesIRs Linear _ _ _ _)) = False

mkMCRGenome :: [Gene] -> [GenesIRsR] -> MultiGIs GenesIRsR
mkMCRGenome = MultiGIs

mkRGenome :: ChromType -> Bool -> Sign -> [Int] -> [Int] -> GenesIRsR
mkRGenome chtype ext sign genes_ irs = GLR $ GenesIRs chtype sign (Vec.fromList . coerce $ genes) (Vec.fromList . map R $ irs) (Gene new)
  where
    genes = if ext then 0 : (genes_ ++ [maxBound]) else genes_
    new = maximum (map abs genes_) + 1

mkRGenome0 :: ChromType -> Sign -> [Int] -> GenesIRsR
mkRGenome0 chtype sign genes = mkRGenome chtype True sign genes (replicate (length genes - 1) 0)

mkFGenome :: ChromType -> Bool -> Sign -> [Int] -> [(Int, Int)] -> GenesIRsF
mkFGenome chtype ext sign genes_ irs = GLF $ GenesIRs chtype sign (Vec.fromList . coerce $ genes) (Vec.fromList . map (uncurry F) $ irs) (Gene new)
  where
    genes = if ext then 0 : (genes_ ++ [maxBound]) else genes_
    new = maximum (map abs genes_) + 1

mkFGenome0 :: ChromType -> Sign -> [Int] -> GenesIRsF
mkFGenome0 chtype sign genes = mkFGenome chtype True sign genes (replicate (length genes - 1) (0, 0))

writeIR :: IR -> BS.ByteString
writeIR ir = LBS.toStrict $ case ir of
  R i -> toLazyByteString . intDec $ i
  F l u -> (toLazyByteString . intDec $ l) <> ":" <> (toLazyByteString . intDec $ u)

splitElements :: BS.ByteString -> [BS.ByteString]
splitElements = filter (not . BS.null) . BS.splitWith (\x -> x == ',' || x == ' ')

readRGenome :: ChromType -> Bool -> Sign -> BS.ByteString -> BS.ByteString -> GenesIRsR
readRGenome chtype extend sign bs_genes bs_irs = mkRGenome chtype extend sign genes irs
  where
    readInt = fst . fromJust . BS.readInt
    genes = map readInt . splitElements $ bs_genes
    irs = map readInt . splitElements $ bs_irs

writeRGenome :: Bool -> GenesIRsR -> (BS.ByteString, BS.ByteString)
writeRGenome rext g@(GLR (GenesIRs _ _ genes irs _)) =
  ( BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec . coerce) $ l_genes,
    BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec) $ l_irs
  )
  where
    l_genes = Vec.toList . (if rext && isLinear g then Vec.slice 1 (size g - 2) else id) $ genes
    l_irs = (\case R i -> i; F _ _ -> error patternError) <$> Vec.toList irs

readFGenome :: ChromType -> Bool -> Sign -> BS.ByteString -> BS.ByteString -> GenesIRsF
readFGenome chtype extend sign bs_genes bs_irs = mkFGenome chtype extend sign genes irs
  where
    readInt = fst . fromJust . BS.readInt
    genes = map readInt . splitElements $ bs_genes
    irs = map readF . splitElements $ bs_irs
    readF fir =
      case BS.splitWith (== ':') fir of
        [lir, uir] -> (readInt lir, readInt uir)
        _ -> error inputError

writeFGenome :: Bool -> GenesIRsF -> (BS.ByteString, BS.ByteString)
writeFGenome rext g@(GLF (GenesIRs _ _ genes irs _)) =
  ( BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec . coerce) $ l_genes,
    BS.unwords . fmap LBS.toStrict $ l_irs
  )
  where
    l_genes = Vec.toList . (if rext && isLinear g then Vec.slice 1 (size g - 2) else id) $ genes
    l_irs = (\case R _ -> error patternError; F l u -> irToBS l u) <$> Vec.toList irs
    irToBS l u = (toLazyByteString . intDec $ l) <> ":" <> (toLazyByteString . intDec $ u)

writeMultiC :: (MultiChromosome g, Chromosome (Chrom g)) => (Chrom g -> (BS.ByteString, BS.ByteString)) -> g -> (BS.ByteString, BS.ByteString)
writeMultiC writer mc = foldr go ("", "") [1 .. mkIdx (numChromosomes mc)]
  where
    go i (bsg, bsi) =
      let chr = getChromosome i mc
          (chr_bsg, chr_bsi) = writer chr
       in (bsg <> (if isLinear chr then "L " else "C ") <> chr_bsg <> " ; ", bsi <> chr_bsi <> " ; ")

flexibilize :: Int -> GenesIRsR -> GenesIRsF
flexibilize l (GLR (GenesIRs chtype sign genes irs new)) = GLF $ GenesIRs chtype sign genes irs' new
  where
    irs' = (\i -> F (i - (l * i `div` 100)) (i + (l * i `div` 100))) . irToInt <$> irs

flexibilizeMC :: Int -> MultiGIs GenesIRsR -> MultiGIs GenesIRsF
flexibilizeMC l (MultiGIs caps_ chrs) = MultiGIs caps_ (map (flexibilize l) chrs)

-- TODO: Adapt operations to allow circular indices
instance RigidIntergenicChromosome GenesIRsR where
  intergenicReversal i j x y (GLR g@(GenesIRs chtype sign vs vi new)) =
    assert (case chtype of { Linear -> 2; Circular -> 1 } <= i)
      . assert (i < j)
      . assert (j <= coerce (case chtype of Linear -> size g - 1; Circular -> size g))
      . assert (0 <= x && x <= ir_x)
      . assert (0 <= y && y <= ir_y)
      . GLR
      $ GenesIRs chtype sign vs' vi' new
    where
      vs' = Vec.modify updateG vs
      vi' = Vec.modify updateIR vi

      updateG :: MVec.MVector s Gene -> ST s ()
      updateG v = do
        mapM_ (\k -> MVec.swap v (coerce $ i + k - 1) (coerce $ j - k - 1)) [0 .. (j - i + 1) `div` 2 - 1]
        case sign of
          Unsigned -> pure ()
          Signed -> mapM_ (MVec.modify v invOri . coerce) [i - 1 .. j - 1]

      updateIR :: MVec.MVector s IR -> ST s ()
      updateIR v = do
        mapM_ (\k -> MVec.swap v (coerce i + coerce k - 1) (coerce j - coerce k - 2)) [0 .. (j - i + 1) `div` 2 - 1]
        MVec.write v (coerce i - 2) (R (x + y))
        MVec.write v (coerce j - 1) (R (x_rest + y_rest))

      ir_x = irToInt (vi ! (coerce i - 2))
      x_rest = ir_x - x
      ir_y = irToInt (vi ! (coerce j - 1))
      y_rest = ir_y - y

  intergenicTransposition i j k x y z (GLR g@(GenesIRs chtype sign vs vi new)) =
    assert (case chtype of { Linear -> 2; Circular -> 1 } <= i)
      . assert (i < j)
      . assert (j < k)
      . assert (k <= coerce (case chtype of Linear -> size g; Circular -> size g + 1))
      . assert (0 <= x && x <= ir_x)
      . assert (0 <= y && y <= ir_y)
      . assert (0 <= z && z <= ir_z)
      . GLR
      $ GenesIRs chtype sign vs' vi' new
    where
      vs' = Vec.modify updateG vs
      vi' = Vec.modify updateIR vi

      updateG :: MVec.MVector s Gene -> ST s ()
      updateG v =
        do
          aux1 <- MVec.clone . MVec.slice (coerce i - 1) (coerce $ j - i) $ v
          aux2 <- MVec.clone . MVec.slice (coerce j - 1) (coerce $ k - j) $ v
          MVec.move (MVec.slice (coerce i - 1) (coerce $ k - j) v) aux2
          MVec.move (MVec.slice (coerce $ i + k - j - 1) (coerce $ j - i) v) aux1

      updateIR :: MVec.MVector s IR -> ST s ()
      updateIR v = do
        do
          aux1 <- MVec.clone . MVec.slice (coerce i - 1) (coerce $ j - i) $ v
          aux2 <- MVec.clone . MVec.slice (coerce j - 1) (coerce $ k - j) $ v
          MVec.move (MVec.slice (coerce i - 1) (coerce $ k - j) v) aux2
          MVec.move (MVec.slice (coerce $ i + k - j - 1) (coerce $ j - i) v) aux1
        MVec.write v (coerce i - 1 - 1) (R (x + y_rest))
        MVec.write v (coerce $ i + k - j - 2) (R (z + x_rest))
        MVec.write v (coerce k - 1 - 1) (R (y + z_rest))

      ir_x = irToInt (vi ! (coerce i - 2))
      x_rest = ir_x - x
      ir_y = irToInt (vi ! (coerce j - 2))
      y_rest = ir_y - y
      ir_z = irToInt (vi ! (coerce k - 2))
      z_rest = ir_z - z

  intergenicDeletion i j x (GLR g@(GenesIRs chtype sign vs vi new)) =
    assert (case chtype of { Linear -> 2; Circular -> 1 } <= i)
      . assert (i < j)
      . assert (j <= coerce (case chtype of Linear -> size g; Circular -> size g + 1))
      . assert (0 <= x && x <= ir_i + ir_j)
      . GLR
      $ GenesIRs chtype sign vs' vi' new
    where
      vs' =
        Vec.slice 0 (coerce i - 1) vs
          Vec.++ Vec.slice (coerce j - 1) (Vec.length vs - coerce j + 1) vs
      vi' =
        Vec.slice 0 (coerce i - 2) vi
          Vec.++ Vec.fromList [R x]
          Vec.++ if coerce j == size g
            then Vec.empty
            else Vec.slice (coerce j - 1) (Vec.length vi - coerce j + 1) vi

      ir_i = irToInt (vi ! (coerce i - 2))
      ir_j = irToInt (vi ! (coerce j - 2))

  -- Genome to be inserted must be open (n genes and n+1 itergenic regions)
  intergenicInsertion i (GLR (GenesIRs chtype _ vsx vix _)) (GLR g@(GenesIRs _ sign vs vi new)) =
    assert (1 <= i)
      . assert (i <= coerce (case chtype of Linear -> size g - 1; Circular -> size g))
      . GLR
      $ GenesIRs chtype sign vs' vi' new'
    where
      vs' =
        Vec.slice 0 (coerce i) vs
          Vec.++ vsx
          Vec.++ Vec.slice (coerce i) (Vec.length vs - coerce i) vs
      new' = Vec.foldr (\gene n -> if abs gene >= n then abs gene + 1 else n) new vsx
      vi' =
        ( if i == 1
            then Vec.empty
            else Vec.slice 0 (coerce i - 1) vi
        )
          Vec.++ vix
          Vec.++ Vec.slice (coerce i) (Vec.length vi - coerce i) vi

instance MultiChromosome (MultiGIs GenesIRsR) where
  type Chrom (MultiGIs GenesIRsR) = GenesIRsR
  getChromosome idx (MultiGIs _ chroms) = chroms !! (idxToInt idx - 1)
  numChromosomes = length . chromosomes

instance MultiChromosome (MultiGIs GenesIRsF) where
  type Chrom (MultiGIs GenesIRsF) = GenesIRsF
  getChromosome idx (MultiGIs _ chroms) = chroms !! (idxToInt idx - 1)
  numChromosomes = length . chromosomes

instance RigidIntergenicMultiChromGenome (MultiGIs GenesIRsR) where
  intergenicDCJ (chri, i) (chrj, j) connectAC x y mc =
    assert (1 <= chri)
      . assert (chri <= chrj)
      . assert (chrj <= mkIdx (numChromosomes mc))
      . assert (1 <= i)
      . assert ((if chri == chrj then i else 0) < j)
      . assert (i <= mkIdx (size chr_i))
      . assert (j <= mkIdx (size chr_j))
      . assert (0 <= x && x <= ir_x)
      . assert (0 <= y && y <= ir_y)
      $ if
          -- ...------i--------...------j--------...
          -- ...------i i+1----...------j j+1----...
          -- ...------i j------ i+1----...j+1----...
          -- ...--------------- i+1----...j+1----...
          -- ...---------------...----i+1 j+1----...
          -- ...---------------...---------------...
          | chri /= chrj && connectAC && isLinear chr_i && isLinear chr_j -> join (chri + 1) (chrj + 1) True True False . join chri chrj' False False True $ mc_cutted
          -- ...(-----i-------)...------j--------...
          -- ...(-----i i+1---)...------j j+1----...
          -- ...i+1-----------i...------j j+1----...
          -- ...i+1-----------i j------...j+1----...
          -- ...i+1--------------------...j+1----...
          -- ...----------------------i+1 j+1----...
          -- ...---------------------------------...
          | chri /= chrj && connectAC && isCircular chr_i && isLinear chr_j -> join chri chrj True True False . join chri chrj' False False True $ mc_cutted
          -- ------i--------...(-----j-------)...
          -- ------i i+1----...(-----j j+1---)...
          -- ------i--------...(-----j-------)...
          -- ------i i+1----...(-----j j+1---)...
          -- ------i i+1----...j+1-----------j...
          -- ------i j-----------j+1 i+1------...
          -- --------------------j+1 i+1------...
          -- ---------------------------------...
          | chri /= chrj && connectAC && isLinear chr_i && isCircular chr_j -> join chri (chri + 1) False False False . join chri chrj' False False True $ mc_cutted
          -- (-----i-------)...(-----j-------)...
          -- (-----i i+1---)...(-----j j+1---)...
          -- i+1-----------i...j+1-----------j...
          -- i+1-----------i j-------------j+1...
          -- i+1---------------------------j+1...
          -- (-------------------------------)...
          | chri /= chrj && connectAC && isCircular chr_i && isCircular chr_j -> turnCircularAtIndex chri . join chri chrj' False False True $ mc_cutted
          -- ------i--------...------j--------...
          -- ------i i+1----...------j j+1----...
          -- ------i j+1---- i+1----...------j...
          -- --------------- i+1----...------j...
          -- --------------- ----i+1 j--------...
          -- --------------- -----------------...
          | chri /= chrj && not connectAC && isLinear chr_i && isLinear chr_j -> join (chri + 1) chrj' False True True . join chri (chrj' + 1) False False False $ mc_cutted
          -- (-----i-------)...------j--------...
          -- (-----i i+1---)...------j j+1----...
          -- i+1-----------i...------j j+1----...
          -- i+1-----------i j+1----...------j...
          -- i+1--------------------...------j...
          -- --------------------i+1 j--------...
          -- ---------------------------------...
          | chri /= chrj && not connectAC && isCircular chr_i && isLinear chr_j -> join chri chrj' False True True . join chri (chrj' + 1) False False False $ mc_cutted
          -- ------i--------...(-----j-------)...
          -- ------i i+1----...(-----j j+1---)...
          -- ------i i+1----...j+1-----------j...
          -- ------i j+1-----------j i+1------...
          -- ----------------------j i+1------...
          -- ---------------------------------...
          | chri /= chrj && not connectAC && isLinear chr_i && isCircular chr_j -> join chri (chri + 1) False False False . join chri chrj' False False False $ mc_cutted
          -- (-----i-------)...(-----j-------)...
          -- (-----i i+1---)...(-----j j+1---)...
          -- i+1-----------i...j+1-----------j...
          -- i+1-----------i j+1-------------j...
          -- i+1-----------------------------j...
          -- (-------------------------------)...
          | chri /= chrj && not connectAC && isCircular chr_i && isCircular chr_j -> turnCircularAtIndex chri . join chri chrj' False False False $ mc_cutted
          -- ------i-----------------j--------...
          -- ------i i+1-------------j j+1----...
          -- ------i j-------------i+1 j+1----...
          -- ----------------------i+1 j+1----...
          -- ---------------------------------...
          | chri == chrj && connectAC && isLinear chr_i -> join chri (chri + 1) False False False . join chri (chri + 1) False False True $ mc_cutted
          -- (-----i-----------------j-------)...
          -- (-----i i+1-------------j-------)...
          -- i+1--------------j--------------i...
          -- i+1--------------j j+1----------i...
          -- j--------------i+1 j+1----------i...
          -- (-------------------------------)...
          | chri == chrj && connectAC && isCircular chr_i -> turnCircularAtIndex chri . join chri (chri + 1) True True False $ mc_cutted
          -- ------i-----------------j--------...
          -- ------i i+1-------------j j+1----...
          -- ------i j+1---- i+1-------------j...
          -- --------------- i+1-------------j...
          -- --------------- (i+1-----------j)...
          | chri == chrj && not connectAC && isLinear chr_i -> turnCircularAtIndex (chri + 1) . join chri (chri + 2) False False False $ mc_cutted
          -- (-----i-----------------j-------)...
          -- (-----i i+1-------------j-------)...
          -- i+1--------------j--------------i...
          -- i+1--------------j j+1----------i...
          -- i+1--------------j (j+1--------i)...
          -- (i+1------------j) (j+1--------i)...
          | chri == chrj && not connectAC && isCircular chr_i -> turnCircularAtIndex chri . turnCircularAtIndex (chri + 1) $ mc_cutted
          | otherwise -> error logicError
    where
      chr_i = getChromosome chri mc
      chr_j = getChromosome chrj mc
      mc_cutted = cut (chrj', j') y . cut (chri, i) x $ mc
      j' = if chri == chrj then j - i + 1 else j
      chrj' = if isLinear chr_i then chrj + 1 else chrj
      ir_x = irToInt (getIR i chr_i)
      ir_y = irToInt (getIR j chr_j)

  cut (chr, i) x mc =
    MultiGIs (caps mc) $
      case splitAt (idxToInt chr - 1) . chromosomes $ mc of
        (_, []) -> error indexError
        (hgs, g : tgs) -> if isLinear g then hgs ++ (g1 : g2 : tgs) else hgs ++ (g' : tgs)
          where
            sig = if hasOri g then Signed else Unsigned
            cap = mkRGenome Linear False sig [geneToInt . head . caps $ mc] []
            ir_idx = getIR i g
            g1 = combine (subGenome 1 i g) (R x) cap
            g2 = combine cap (sumIR ir_idx (R (-x))) (subGenome (i + 1) (mkIdx (size g)) g)
            g' = combine cap (sumIR ir_idx (R (-x))) (combine (turnLinear i g) (R x) cap)
  join chri chrj move_first invi invj mc =
    if chri > chrj
      then join chrj chri move_first invj invi mc
      else
        MultiGIs (caps mc) $
          let (gs1, gs23) = splitAt (idxToInt chri - 1) . chromosomes $ mc
           in case splitAt (idxToInt (chrj - chri)) gs23 of
                (_, []) -> error indexError
                ([], _) -> error indexError
                (gi : gs2, gj : gs3) -> gs1 ++ (if move_first then gs2 ++ g' : gs3 else g' : gs2 ++ gs3)
                  where
                    g' = combine (subGenome 1 (mkIdx $ size gi' - 1) gi') (R $ irToInt (getIR (mkIdx $ size gi' - 1) gi') + irToInt (getIR 1 gj')) (subGenome 2 (mkIdx $ size gj') gj')
                    gi' = if invi then intergenicFullReversal gi else gi
                    gj' = if invj then intergenicFullReversal gj else gj

  turnCircularAtIndex chri mc = MultiGIs (caps mc) . replace (idxToInt chri - 1) turnCircular . chromosomes $ mc

class Matcher m g1 g2 where
  isMatch :: m g1 g2 -> g1 -> g2 -> Bool
  isDirectMatch :: m g1 g2 -> g1 -> g2 -> Bool
  isReverseMatch :: m g1 g2 -> g1 -> g2 -> Bool
  areBalanced :: m g1 g2 -> g1 -> g2 -> Bool

newtype FlipMatcher m g1 g2 = FlipMatcher m

data RigidRigidDirectMatcher g1 g2 = RRDM

data RigidFlexibleDirectMatcher g1 g2 = RFDM

data RigidRigidReverseMatcher g1 g2 = RRRM

data RigidFlexibleReverseMatcher g1 g2 = RFRM

instance (Matcher m g1 g2) => Matcher (FlipMatcher (m g1 g2)) g2 g1 where
  isMatch (FlipMatcher matcher) g1 g2 = isMatch matcher g2 g1
  isDirectMatch (FlipMatcher matcher) g1 g2 = isDirectMatch matcher g2 g1
  isReverseMatch (FlipMatcher matcher) g1 g2 = isReverseMatch matcher g2 g1
  areBalanced (FlipMatcher matcher) g1 g2 = areBalanced matcher g2 g1

instance Matcher RigidRigidDirectMatcher GenesIRsR GenesIRsR where
  isMatch = isDirectMatch
  isDirectMatch _ g h = g == h
  isReverseMatch _ _ _ = False
  areBalanced _ (GLR (GenesIRs _ _ genesG irsG _)) (GLR (GenesIRs _ _ genesH irsH _)) =
    List.sort (fmap canonicOri (Vec.toList genesG)) == List.sort (fmap canonicOri (Vec.toList genesH))
      && sum (fmap irToInt irsG) == sum (fmap irToInt irsH)

instance Matcher RigidFlexibleDirectMatcher GenesIRsR GenesIRsF where
  isMatch = isDirectMatch
  isDirectMatch RFDM (GLR (GenesIRs Linear _ genesG irsG _)) (GLF (GenesIRs Linear _ genesH irsH _)) =
    genesG == genesH && and (Vec.zipWith check irsG irsH)
    where
      check :: IR -> IR -> Bool
      check (R ir) (F irl irr) = (irl <= ir) && (ir <= irr)
      check _ _ = error patternError
  isDirectMatch RFDM (GLR (GenesIRs Circular _ _ _ _)) (GLF _) = undefined
  isDirectMatch RFDM (GLR _) (GLF (GenesIRs Circular _ _ _ _)) = undefined
  isReverseMatch _ _ _ = False
  areBalanced _ (GLR (GenesIRs _ _ genesG irsG _)) (GLF (GenesIRs _ _ genesH irsH _)) =
    List.sort (fmap canonicOri (Vec.toList genesG)) == List.sort (fmap canonicOri (Vec.toList genesH))
      && sum (fmap irToIntLow irsH) <= sum (fmap irToInt irsG)
      && sum (fmap irToInt irsG) <= sum (fmap irToIntHigh irsH)
    where
      irToIntLow ir = case ir of
        R _ -> error patternError
        F l _ -> l
      irToIntHigh ir = case ir of
        R _ -> error patternError
        F _ r -> r

instance Matcher RigidRigidReverseMatcher GenesIRsR GenesIRsR where
  isMatch RRRM g h = isDirectMatch RRRM g h || isReverseMatch RRRM g h
  isDirectMatch _ = isDirectMatch RRDM
  isReverseMatch _ g = isDirectMatch RRDM (intergenicFullReversal g)
  areBalanced _ = areBalanced RRDM

instance Matcher RigidFlexibleReverseMatcher GenesIRsR GenesIRsF where
  isMatch RFRM g h = isDirectMatch RFRM g h || isReverseMatch RFRM g h
  isDirectMatch _ = isDirectMatch RFDM
  isReverseMatch _ g = isDirectMatch RFDM (intergenicFullReversal g)
  areBalanced _ = areBalanced RFDM

-- | check if @g1@ is compatible with a subgenome of @g2@ and return the index
-- of the first gene from the first compatible subgenome
isCompatibleWithSubGenome :: (Matcher m g1 g2, Genome g1, Genome g2) => m g1 g2 -> g1 -> g2 -> Maybe Idx
isCompatibleWithSubGenome matcher g1 g2 = find testSubgenome [1 .. mkIdx (size g2)]
  where
    testSubgenome beg =
      let end = beg + mkIdx (size g1) - 1
       in end <= mkIdx (size g2) && isMatch matcher g1 (subGenome beg end g2)

randomGenome :: MonadRandom mon => Bool -> Int -> Int -> Sign -> mon GenesIRsR
randomGenome zeros n lim signed = do
  coins <- getRandoms
  ls <- case signed of
    Unsigned -> take n <$> getRandomRs (1 :: Gene, coerce lim)
    Signed -> zipWith swaps coins . take n <$> getRandomRs (1 :: Gene, coerce lim)
  li <- take (n + 1) <$> if zeros then return (repeat 0) else getRandomRs (1, 100)
  return $ mkRGenome Linear True signed (coerce ls) li
  where
    swaps b v = if b then v else invOri v

randomGenomeWithReplicas :: MonadRandom mon => Bool -> Int -> Int -> Int -> Int -> Sign -> mon GenesIRsR
randomGenomeWithReplicas zeros n rep low high signed = do
  coins <- getRandoms
  occs <- getRandomRs (low, high)
  let ls =
        ( case signed of
            Unsigned -> id
            Signed -> zipWith swaps coins
        )
          . (\l -> l ++ coerce [length l + 1 .. n])
          . concat
          . zipWith replicate occs
          $ [Gene 1 .. Gene rep]
  li <- take (n + 1) <$> if zeros then return (repeat 0) else getRandomRs (1, 100)
  return $ mkRGenome Linear True signed (coerce ls) li
  where
    swaps b v = if b then v else invOri v

shuffleGenome :: (MonadRandom mon) => GenesIRsR -> mon GenesIRsR
shuffleGenome g@(GLR (GenesIRs chtype sign ls_ li _)) = do
  coins <- getRandoms
  let s = sum . fmap irToInt $ li
      ls = Vec.toList ls_
  ls' <- case sign of
    Unsigned -> shuffleM ls
    Signed -> zipWith swaps coins <$> shuffleM ls
  x <- List.sort . take (size g) <$> getRandomRs (0, s)
  let li' = zipWith (-) (x ++ [s]) (0 : x)
  return $ mkRGenome chtype True sign (coerce ls') li'
  where
    swaps b v = if b then v else invOri v
