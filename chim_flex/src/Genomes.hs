{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE ExplicitForAll #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE OverloadedStrings #-}

module Genomes
  ( Genome (..),
    IntergenicGenome (..),
    RigidIntergenicGenome (..),
    Gene,
    mkGene,
    geneToInt,
    IR,
    irToInt,
    Idx,
    Sign (..),
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
    GenesIRsR,
    GenesIRsF,
    mkRGenome,
    mkFGenome,
    mkRGenome0,
    mkFGenome0,
    writeIR,
    readRGenome,
    writeRGenome,
    readFGenome,
    writeFGenome,
    flexibilize,
    Matcher (..),
    RigidRigidDirectMatcher (..),
    RigidFlexibleDirectMatcher (..),
    RigidRigidReverseMatcher (..),
    RigidFlexibleReverseMatcher (..),
    randomGenome,
    randomGenomeWithReplicas,
    shuffleGenome,
  )
where

import Control.Exception (assert)
import Control.Monad.Random (MonadRandom, getRandomRs, getRandoms)
import Data.ByteString.Builder (intDec, toLazyByteString)
import Data.ByteString.Char8 qualified as BS
import Data.ByteString.Lazy qualified as LBS
import Data.Coerce (coerce)
import Data.Hashable (Hashable)
import Data.IntMap (IntMap)
import Data.IntMap qualified as IntMap
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

data IR = R Int | F Int Int deriving (Eq)

irToInt :: IR -> Int
irToInt ir = case ir of
  R i -> i
  F _ _ -> error patternError

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

data Sign = Signed | Unsigned deriving (Eq, Show, Enum)

-- Representation with a list of genes and a list of intergenic regions
-- (only ints for regid, or intervals for flexible)
-- We assume that the genome starts and end with genes
data GenesIRs = GenesIRs Sign (Vector Gene) (Vector IR)

instance Show GenesIRs where
  show (GenesIRs _ genes irs) =
    unwords . (("(" ++ head str_s ++ ")") :) $
      zipWith (\ir a -> "- " ++ ir ++ " - (" ++ a ++ ")") str_i (tail str_s)
    where
      str_s = Vec.toList $ (\i -> if i == maxBound then "inf" else show i) <$> genes
      str_i = Vec.toList $ show <$> irs

class Genome g where
  isGene :: g -> Gene -> Bool
  size :: g -> Int

  -- get gene at position i (index starts in 1), requires 1 <= i <= size g
  getGene :: Idx -> g -> Gene

  -- for indices i j, it gets the subgenome starting at gene position i
  -- and ending at gene position j
  -- requires 1 <= i < j <= size g
  subGenome :: Idx -> Idx -> g -> g

  -- Return set with all the gene values
  alphabet :: g -> Set Gene

class (Genome g) => IntergenicGenome g where
  getIR :: Idx -> g -> IR

class (Genome g) => RigidIntergenicGenome g where
  intergenicFullReversal :: g -> g
  intergenicReversal :: Idx -> Idx -> Int -> Int -> g -> g
  intergenicTransposition :: Idx -> Idx -> Idx -> Int -> Int -> Int -> g -> g
  intergenicInsertion :: Idx -> g -> g -> g
  intergenicDeletion :: Idx -> Idx -> Int -> g -> g

occurrence :: (Genome g) => g -> Gene -> Int
occurrence g a = sum . fmap (\i -> if abs (getGene (Idx i) g) == a then 1 else 0) $ [1 .. size g]

newtype GeneMap a = GM (IntMap a) deriving newtype (Functor, Foldable)

positionMap :: (Genome g) => g -> GeneMap [Idx]
positionMap g = GM . IntMap.fromListWith (++) . map (\i -> (geneToInt . abs $ getGene i g, [i])) $ [1 .. mkIdx (size g)]

geneMapLookup :: Gene -> GeneMap a -> Maybe a
geneMapLookup a (GM gm) = IntMap.lookup (geneToInt $ abs a) gm

geneMapAdjust :: (a -> a) -> Gene -> GeneMap a -> GeneMap a
geneMapAdjust f a (GM gm) = GM $ IntMap.adjust f (geneToInt $ abs a) gm

occurrenceMax :: (Genome g) => g -> Int
occurrenceMax = maximum . fmap length . positionMap

instance Genome GenesIRs where
  isGene (GenesIRs _ genes _) a = a `elem` genes
  size (GenesIRs _ genes _) = Vec.length genes
  getGene idx (GenesIRs _ genes _) = genes Vec.! (idxToInt idx - 1)
  subGenome idxi idxj (GenesIRs sign genes irs) = GenesIRs sign (Vec.slice i n genes) (Vec.slice i (n - 1) irs)
    where
      i = idxToInt idxi - 1
      j = idxToInt idxj
      n = j - i
  alphabet (GenesIRs _ genes _) = Set.fromList . Vec.toList . fmap abs $ genes

instance IntergenicGenome GenesIRs where
  getIR idx (GenesIRs _ _ irs) = irs Vec.! (idxToInt idx - 1)

newtype GenesIRsR = GLR GenesIRs deriving newtype (Show, Genome, IntergenicGenome)

newtype GenesIRsF = GLF GenesIRs deriving newtype (Show, Genome, IntergenicGenome)

instance Eq GenesIRsR where
  (GLR (GenesIRs sign1 genes1 irs1)) == (GLR (GenesIRs sign2 genes2 irs2)) =
    sign1 == sign2 && genes1 == genes2 && irs1 == irs2

mkRGenome :: Bool -> Sign -> [Int] -> [Int] -> GenesIRsR
mkRGenome ext sign genes_ irs = GLR $ GenesIRs sign (Vec.fromList . coerce $ genes) (Vec.fromList . map R $ irs)
  where
    genes = if ext then 0 : (genes_ ++ [maxBound]) else genes_

mkRGenome0 :: Sign -> [Int] -> GenesIRsR
mkRGenome0 sign genes = mkRGenome True sign genes (replicate (length genes - 1) 0)

mkFGenome :: Bool -> Sign -> [Int] -> [(Int, Int)] -> GenesIRsF
mkFGenome ext sign genes_ irs = GLF $ GenesIRs sign (Vec.fromList . coerce $ genes) (Vec.fromList . map (uncurry F) $ irs)
  where
    genes = if ext then 0 : (genes_ ++ [maxBound]) else genes_

mkFGenome0 :: Sign -> [Int] -> GenesIRsF
mkFGenome0 sign genes = mkFGenome True sign genes (replicate (length genes - 1) (0, 0))

writeIR :: IR -> BS.ByteString
writeIR ir = LBS.toStrict $ case ir of
  R i -> toLazyByteString . intDec $ i
  F l u -> (toLazyByteString . intDec $ l) <> ":" <> (toLazyByteString . intDec $ u)

readRGenome :: Bool -> Sign -> BS.ByteString -> BS.ByteString -> GenesIRsR
readRGenome extend sign bs_genes bs_irs = mkRGenome extend sign genes irs
  where
    readInt = fst . fromJust . BS.readInt
    genes = map readInt . BS.splitWith (\x -> x == ',' || x == ' ') $ bs_genes
    irs = map readInt . BS.splitWith (\x -> x == ',' || x == ' ') $ bs_irs

writeRGenome :: Bool -> GenesIRsR -> (BS.ByteString, BS.ByteString)
writeRGenome rext g@(GLR (GenesIRs _ genes irs)) =
  ( BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec . coerce) $ l_genes,
    BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec) $ l_irs
  )
  where
    l_genes = Vec.toList . (if rext then Vec.slice 1 (size g - 2) else id) $ genes
    l_irs = (\case R i -> i; F _ _ -> error patternError) <$> Vec.toList irs

readFGenome :: Bool -> Sign -> BS.ByteString -> BS.ByteString -> GenesIRsF
readFGenome extend sign bs_genes bs_irs = mkFGenome extend sign genes irs
  where
    readInt = fst . fromJust . BS.readInt
    genes = map readInt . BS.splitWith (\x -> x == ',' || x == ' ') $ bs_genes
    irs = map readF . BS.splitWith (\x -> x == ',' || x == ' ') $ bs_irs
    readF fir =
      case BS.splitWith (== ':') fir of
        [lir, uir] -> (readInt lir, readInt uir)
        _ -> error patternError

writeFGenome :: Bool -> GenesIRsF -> (BS.ByteString, BS.ByteString)
writeFGenome rext g@(GLF (GenesIRs _ genes irs)) =
  ( BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec . coerce) $ l_genes,
    BS.unwords . fmap LBS.toStrict $ l_irs
  )
  where
    l_genes = Vec.toList . (if rext then Vec.slice 1 (size g - 2) else id) $ genes
    l_irs = (\case R _ -> error patternError; F l u -> irToBS l u) <$> Vec.toList irs
    irToBS l u = (toLazyByteString . intDec $ l) <> ":" <> (toLazyByteString . intDec $ u)

flexibilize :: Int -> GenesIRsR -> GenesIRsF
flexibilize l (GLR (GenesIRs sign genes irs)) = GLF $ GenesIRs sign genes irs'
  where
    irs' = (\i -> F (i - (l * i `div` 100)) (i + (l * i `div` 100))) . irToInt <$> irs

instance RigidIntergenicGenome GenesIRsR where
  intergenicFullReversal (GLR (GenesIRs sign vs vi)) = GLR $ GenesIRs sign vs' vi'
    where
      vs' = case sign of { Unsigned -> id; Signed -> fmap negate } $ Vec.reverse vs
      vi' = Vec.reverse vi

  intergenicReversal i j x y (GLR g@(GenesIRs sign vs vi)) =
    assert (2 <= i)
      . assert (i < j)
      . assert (j <= coerce (size g) - 1)
      . assert (0 <= x && x <= ir_x)
      . assert (0 <= y && y <= ir_y)
      . GLR
      $ GenesIRs sign vs' vi'
    where
      vs' = Vec.modify updateG vs
      vi' = Vec.modify updateIR vi

      updateG v = do
        mapM_ (\k -> MVec.swap v (coerce $ i + k - 1) (coerce $ j - k - 1)) [0 .. (j - i + 1) `div` 2 - 1]
        case sign of
          Unsigned -> pure ()
          Signed -> mapM_ (MVec.modify v invOri . coerce) [i .. j]
        mapM_ (MVec.modify v invOri . coerce) [i .. j]

      updateIR v = do
        mapM_ (\k -> MVec.swap v (coerce i + coerce k - 1) (coerce j - coerce k - 2)) [0 .. (j - i + 1) `div` 2 - 1]
        MVec.write v (coerce i - 2) (R (x + y))
        MVec.write v (coerce j - 1) (R (x_rest + y_rest))

      ir_x = irToInt (vi ! (coerce i - 2))
      x_rest = ir_x - x
      ir_y = irToInt (vi ! (coerce j - 1))
      y_rest = ir_y - y

  intergenicTransposition i j k x y z (GLR g@(GenesIRs sign vs vi)) =
    assert (2 <= i)
      . assert (i < j)
      . assert (j < k)
      . assert (k <= Idx (size g))
      . assert (0 <= x && x <= ir_x)
      . assert (0 <= y && y <= ir_y)
      . assert (0 <= z && z <= ir_z)
      . GLR
      $ GenesIRs sign vs' vi'
    where
      vs' = Vec.modify updateG vs
      vi' = Vec.modify updateIR vi

      updateG v =
        do
          aux1 <- MVec.clone . MVec.slice (coerce i - 1) (coerce $ j - i) $ v
          aux2 <- MVec.clone . MVec.slice (coerce j - 1) (coerce $ k - j) $ v
          MVec.move (MVec.slice (coerce i - 1) (coerce $ k - j) v) aux2
          MVec.move (MVec.slice (coerce $ i + k - j - 1) (coerce $ j - i) v) aux1

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

  intergenicDeletion i j x (GLR g@(GenesIRs sign vs vi)) =
    assert (2 <= i)
      . assert (i < j)
      . assert (j <= coerce (size g))
      . assert (0 <= x && x <= ir_i + ir_j)
      . GLR
      $ GenesIRs sign vs' vi'
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
  intergenicInsertion i (GLR (GenesIRs _ vsx vix)) (GLR g@(GenesIRs sign vs vi)) =
    assert (1 <= i)
      . assert (i <= coerce (size g - 1))
      . GLR
      $ GenesIRs sign vs' vi'
    where
      vs' =
        Vec.slice 0 (coerce i) vs
          Vec.++ vsx
          Vec.++ Vec.slice (coerce i) (Vec.length vs - coerce i) vs
      vi' =
        ( if i == 1
            then Vec.empty
            else Vec.slice 0 (coerce i - 1) vi
        )
          Vec.++ vix
          Vec.++ Vec.slice (coerce i) (Vec.length vi - coerce i) vi

class Matcher m g1 g2 where
  isMatch :: m g1 g2 -> g1 -> g2 -> Bool
  isDirectMatch :: m g1 g2 -> g1 -> g2 -> Bool
  isReverseMatch :: m g1 g2 -> g1 -> g2 -> Bool

data RigidRigidDirectMatcher g1 g2 = RRDM

data RigidFlexibleDirectMatcher g1 g2 = RFDM

data RigidRigidReverseMatcher g1 g2 = RRRM

data RigidFlexibleReverseMatcher g1 g2 = RFRM

instance Matcher RigidRigidDirectMatcher GenesIRsR GenesIRsR where
  isMatch = isDirectMatch
  isDirectMatch _ g h = g == h
  isReverseMatch _ g = isMatch RRDM (intergenicFullReversal g)

instance Matcher RigidFlexibleDirectMatcher GenesIRsR GenesIRsF where
  isMatch = isDirectMatch
  isDirectMatch RFDM (GLR (GenesIRs _ genesG irsG)) (GLF (GenesIRs _ genesH irsH)) =
    genesG == genesH && and (Vec.zipWith check irsG irsH)
    where
      check :: IR -> IR -> Bool
      check (R ir) (F irl irr) = (irl <= ir) && (ir <= irr)
      check _ _ = error patternError
  isReverseMatch _ g = isMatch RFDM (intergenicFullReversal g)

instance Matcher RigidRigidReverseMatcher GenesIRsR GenesIRsR where
  isMatch RRRM g h = isDirectMatch RRRM g h || isReverseMatch RRRM g h
  isDirectMatch _ = isDirectMatch RRDM
  isReverseMatch _ = isReverseMatch RRDM

instance Matcher RigidFlexibleReverseMatcher GenesIRsR GenesIRsF where
  isMatch RFRM g h = isDirectMatch RFRM g h || isReverseMatch RFRM g h
  isDirectMatch _ = isDirectMatch RFDM
  isReverseMatch _ = isReverseMatch RFDM

randomGenome :: MonadRandom mon => Bool -> Int -> Int -> Sign -> mon GenesIRsR
randomGenome zeros n lim signed = do
  coins <- getRandoms
  ls <- case signed of
    Unsigned -> take n <$> getRandomRs (1 :: Gene, coerce lim)
    Signed -> zipWith swaps coins . take n <$> getRandomRs (1 :: Gene, coerce lim)
  li <- take (n + 1) <$> if zeros then return (repeat 0) else getRandomRs (1, 100)
  return $ mkRGenome True signed (coerce ls) li
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
  return $ mkRGenome True signed (coerce ls) li
  where
    swaps b v = if b then v else invOri v

shuffleGenome :: (MonadRandom mon) => GenesIRsR -> mon GenesIRsR
shuffleGenome g@(GLR (GenesIRs sign ls_ li)) = do
  coins <- getRandoms
  let s = sum . fmap irToInt $ li
      ls = Vec.toList ls_
  ls' <- case sign of
    Unsigned -> shuffleM ls
    Signed -> zipWith swaps coins <$> shuffleM ls
  x <- List.sort . take (size g) <$> getRandomRs (0, s)
  let li' = zipWith (-) (x ++ [s]) (0 : x)
  return $ mkRGenome True sign (coerce ls') li'
  where
    swaps b v = if b then v else invOri v
