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
    IR,
    Idx,
    mkIdx,
    idxDist,
    idxToInt,
    incIdx,
    decIdx,
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
    Matcher (..),
    RigidRigidMatcher (..),
    RigidFlexibleMatcher (..),
  )
where

import Control.Exception (assert)
import Data.ByteString.Builder (intDec, toLazyByteString)
import Data.ByteString.Char8 qualified as BS
import Data.ByteString.Lazy qualified as LBS
import Data.Coerce (coerce)
import Data.Hashable (Hashable)
import Data.Maybe (fromJust)
import Data.Vector (Vector, (!))
import Data.Vector qualified as Vec
import Data.Vector.Mutable qualified as MVec
import LocalBase

newtype Gene = Gene Int deriving newtype (Eq, Show, Read, Hashable, Ord, Num, Bounded, Enum)

mkGene :: Int -> Gene
mkGene = Gene

data IR = R Int | F Int Int deriving (Eq)

newtype Idx = Idx Int deriving newtype (Eq, Show, Read, Hashable, Ord, Num, Bounded, Enum, Integral, Real)

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

-- Representation with a list of genes and a list of intergenic regions
-- (only ints for regid, or intervals for flexible)
-- We assume that the genome starts and end with genes
data GenesIRs = GenesIRs (Vector Gene) (Vector IR)

instance Show GenesIRs where
  show (GenesIRs genes irs) =
    unwords . (("(" ++ head str_s ++ ")") :) $
      zipWith (\ir a -> "- " ++ ir ++ " - (" ++ a ++ ")") str_i (tail str_s)
    where
      str_s = Vec.toList $ (\i -> if i == maxBound then "inf" else show i) <$> genes
      str_i = Vec.toList $ show <$> irs

class (Show g) => Genome g where
  isGene :: g -> Gene -> Bool
  size :: g -> Int

  -- get gene at position i (index starts in 1), requires 1 <= i <= size g
  getGene :: Idx -> g -> Gene

  -- for indices i j, it gets the subgenome starting at gene position i
  -- and ending at gene position j
  -- requires 1 <= i < j <= size g
  subGenome :: Idx -> Idx -> g -> g

class (Genome g) => IntergenicGenome g where
  getIR :: Idx -> g -> IR

class (Genome g) => RigidIntergenicGenome g where
  intergenicReversal :: Idx -> Idx -> Int -> Int -> g -> g

instance Genome GenesIRs where
  isGene (GenesIRs genes _) a = a `elem` genes
  size (GenesIRs genes _) = Vec.length genes
  getGene idx (GenesIRs genes _) = genes Vec.! (idxToInt idx - 1)
  subGenome idxi idxj (GenesIRs genes irs) = GenesIRs (Vec.slice i n genes) (Vec.slice i (n - 1) irs)
    where
      i = idxToInt idxi - 1
      j = idxToInt idxj
      n = j - i

instance IntergenicGenome GenesIRs where
  getIR idx (GenesIRs _ irs) = irs Vec.! (idxToInt idx - 1)

newtype GenesIRsR = GLR GenesIRs deriving newtype (Show, Genome, IntergenicGenome)

newtype GenesIRsF = GLF GenesIRs deriving newtype (Show, Genome, IntergenicGenome)

mkRGenome :: [Int] -> [Int] -> GenesIRsR
mkRGenome genes irs = GLR $ GenesIRs (Vec.fromList . coerce $ genes) (Vec.fromList . map R $ irs)

mkRGenome0 :: [Int] -> GenesIRsR
mkRGenome0 genes = mkRGenome genes (replicate (length genes - 1) 0)

mkFGenome :: [Int] -> [(Int, Int)] -> GenesIRsF
mkFGenome genes irs = GLF $ GenesIRs (Vec.fromList . coerce $ genes) (Vec.fromList . map (uncurry F) $ irs)

mkFGenome0 :: [Int] -> GenesIRsF
mkFGenome0 genes = mkFGenome genes (replicate (length genes - 1) (0, 0))

writeIR :: IR -> BS.ByteString
writeIR ir = LBS.toStrict $ case ir of
  R i -> toLazyByteString . intDec $ i
  F l u -> (toLazyByteString . intDec $ l) <> ":" <> (toLazyByteString . intDec $ u)

readRGenome :: Bool -> BS.ByteString -> BS.ByteString -> GenesIRsR
readRGenome extend bs_genes bs_irs = mkRGenome genes irs
  where
    readInt = fst . fromJust . BS.readInt
    genes_ = map readInt . BS.splitWith (\x -> x == ',' || x == ' ') $ bs_genes
    genes = if extend then 0 : (genes_ ++ [maxBound]) else genes_
    irs = map readInt . BS.splitWith (\x -> x == ',' || x == ' ') $ bs_irs

writeRGenome :: Bool -> GenesIRsR -> (BS.ByteString, BS.ByteString)
writeRGenome rext g@(GLR (GenesIRs genes irs)) =
  ( BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec . coerce) $ l_genes,
    BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec) $ l_irs
  )
  where
    l_genes = Vec.toList . (if rext then Vec.slice 1 (size g - 2) else id) $ genes
    l_irs = (\case R i -> i; F _ _ -> error patternError) <$> Vec.toList irs

readFGenome :: Bool -> BS.ByteString -> BS.ByteString -> GenesIRsF
readFGenome extend bs_genes bs_irs = mkFGenome genes irs
  where
    readInt = fst . fromJust . BS.readInt
    genes_ = map readInt . BS.splitWith (\x -> x == ',' || x == ' ') $ bs_genes
    genes = if extend then 0 : (genes_ ++ [maxBound]) else genes_
    irs = map readF . BS.splitWith (\x -> x == ',' || x == ' ') $ bs_irs
    readF fir =
      case BS.splitWith (== ':') fir of
        [lir, uir] -> (readInt lir, readInt uir)
        _ -> error patternError

writeFGenome :: Bool -> GenesIRsF -> (BS.ByteString, BS.ByteString)
writeFGenome rext g@(GLF (GenesIRs genes irs)) =
  ( BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec . coerce) $ l_genes,
    BS.unwords . fmap LBS.toStrict $ l_irs
  )
  where
    l_genes = Vec.toList . (if rext then Vec.slice 1 (size g - 2) else id) $ genes
    l_irs = (\case R _ -> error patternError; F l u -> irToBS l u) <$> Vec.toList irs
    irToBS l u = (toLazyByteString . intDec $ l) <> ":" <> (toLazyByteString . intDec $ u)

instance RigidIntergenicGenome GenesIRsR where
  intergenicReversal i j x y (GLR g@(GenesIRs vs vi)) =
    assert (2 <= i)
      . assert (i < j)
      . assert (j <= coerce (size g) - 1)
      . assert (0 <= x && x <= ir_x)
      . assert (0 <= y && y <= ir_y)
      . GLR
      $ GenesIRs vs' vi'
    where
      vs' = Vec.modify updateG vs
      vi' = Vec.modify updateIR vi

      updateG v = do
        mapM_ (\k -> MVec.swap v (coerce $ i + k - 1) (coerce $ j - k - 1)) [0 .. (j - i + 1) `div` 2 - 1]
        -- case genomeIsSigned g of
        --   Unsigned -> pure ()
        --   Signed -> mapM_ (MVec.modify v invOri . coerce) [i .. j]
        mapM_ (MVec.modify v invOri . coerce) [i .. j]

      updateIR v = do
        mapM_ (\k -> MVec.swap v (coerce i + coerce k - 1) (coerce j - coerce k - 2)) [0 .. (j - i + 1) `div` 2 - 1]
        MVec.write v (coerce i - 2) (R (x + y))
        MVec.write v (coerce j - 1) (R (x_rest + y_rest))

      ir_x = case vi ! (coerce i - 2) of
        R ir -> ir
        F _ _ -> error patternError
      x_rest = ir_x - x
      ir_y = case vi ! (coerce j - 1) of
        R ir -> ir
        F _ _ -> error patternError
      y_rest = ir_y - y

class Matcher m g1 g2 where
  isMatch :: m g1 g2 -> g1 -> g2 -> Bool

data RigidRigidMatcher g1 g2 = RRM

data RigidFlexibleMatcher g1 g2 = RFM

instance Matcher RigidRigidMatcher GenesIRsR GenesIRsR where
  isMatch RRM (GLR (GenesIRs genesG irsG)) (GLR (GenesIRs genesH irsH)) =
    genesG == genesH && irsG == irsH

instance Matcher RigidFlexibleMatcher GenesIRsR GenesIRsF where
  isMatch RFM (GLR (GenesIRs genesG irsG)) (GLF (GenesIRs genesH irsH)) =
    genesG == genesH && and (Vec.zipWith check irsG irsH)
    where
      check :: IR -> IR -> Bool
      check (R ir) (F irl irr) = (irl <= ir) && (ir <= irr)
      check _ _ = error patternError