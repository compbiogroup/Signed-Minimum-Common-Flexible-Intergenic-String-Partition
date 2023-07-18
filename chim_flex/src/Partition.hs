{-# LANGUAGE ExplicitForAll #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}

module Partition
  ( writePartition,
    getBlocksCorrespondence,
    trivialPartition,
    Partition,
    Breakpoints,
    mkPartition,
    underlineGenome,
    breakpoints,
    blocks,
    checkCommon,
    CommonPartition,
    mkCommonPartition,
    mkCommonPartition2,
    combine,
    blockDelsToBps,
    bpsToBlockDels,
  )
where

import Control.Exception (assert)
import Data.ByteString qualified as BS
import Data.ByteString.Char8 qualified as BS
import Data.EnumSet (EnumSet)
import Data.EnumSet qualified as EnumSet
import Data.Foldable (foldrM, toList)
import Data.Map qualified as Map
import Data.Maybe (isJust)
import Data.Sequence (Seq (Empty, (:<|), (:|>)), (><))
import Data.Sequence qualified as Seq
import Data.Set qualified as Set
import Genomes (GenesIRsF, GenesIRsR, Genome (..), IR, Idx, IntergenicGenome (..), Matcher (..), geneMapLookup, incIdx, mkIdx, positionMap, writeFGenome, writeIR, writeRGenome)
import LocalBase

type Breakpoints = EnumSet Idx

data Partition g where
  GP :: (Genome g) => g -> Breakpoints -> Partition g

mkPartition :: (Genome g) => g -> Breakpoints -> Partition g
mkPartition = GP

underlineGenome :: Partition p -> p
underlineGenome (GP g _) = g

trivialPartition :: Genome g => g -> Partition g
trivialPartition g = GP g EnumSet.empty

breakpoints :: Partition g -> Breakpoints
breakpoints (GP _ bps) = bps

breakpointsIdx :: Partition g -> [Idx]
breakpointsIdx (GP _ bps) = EnumSet.toAscList bps

breakpointsIR :: (IntergenicGenome g) => Partition g -> [IR]
breakpointsIR gp@(GP g _) = map (`getIR` g) (breakpointsIdx gp)

blocks :: Partition g -> [g]
blocks pg@(GP g _) = zipWith (\i succi -> subGenome (incIdx i) succi g) bg $ tail bg
  where
    bg = 0 : breakpointsIdx pg ++ [mkIdx (size g)]

data CommonPartition g1 g2 where
  CGP :: (Genome g1, Genome g2) => Partition g1 -> Partition g2 -> CommonPartition g1 g2

mkCommonPartition :: (Matcher m g1 g2, Genome g1, Genome g2) => m g1 g2 -> Partition g1 -> Partition g2 -> CommonPartition g1 g2
mkCommonPartition matcher pg ph = assert (checkCommon_ matcher pg ph) $ CGP pg ph

-- Verify if the the elements of pg can be assign to elements of ph with a perfect
-- match.
checkCommon :: (Matcher m g1 g2) => m g1 g2 -> CommonPartition g1 g2 -> Bool
checkCommon matcher (CGP pg ph) = checkCommon_ matcher pg ph

checkCommon_ :: (Matcher m g1 g2) => m g1 g2 -> Partition g1 -> Partition g2 -> Bool
checkCommon_ matcher pg ph = isJust $ testPerfectMatch cors [0 .. length cors - 1]
  where
    cors = getBlocksCorrespondence_ matcher pg ph
    testPerfectMatch aux_corres = foldrM (selectBlock Set.empty) Map.empty
      where
        -- Select block from G that corresponds to block represented by v.
        -- fixCorres saves a map between blocks of G and their fix correspondent
        -- in H.
        selectBlock seen v fixCorres =
          if done then Just fixCorres' else Nothing
          where
            (_, fixCorres', done) = selectBlock' seen v fixCorres

        selectBlock' seen v = testCandidates v cands seen
          where
            cands = aux_corres !! v

        testCandidates _ [] seen fixCorres = (seen, fixCorres, False)
        testCandidates v (u : us) seen fixCorres =
          if
              | u `Set.member` seen -> testCandidates v us seen fixCorres
              | u `Map.notMember` fixCorres -> (seen', fixCorres', True)
              | done -> (seen_rec, blocksFromG_rec', True)
              | otherwise -> testCandidates v us seen_rec blocksFromG_rec
          where
            seen' = Set.insert u seen
            fixCorres' = Map.insert u v fixCorres
            (seen_rec, blocksFromG_rec, done) = selectBlock' seen' (fixCorres Map.! u) fixCorres
            blocksFromG_rec' = Map.insert u v blocksFromG_rec

mkCommonPartition2 :: (Matcher m g1 g2, Genome g1, Genome g2) => m g1 g2 -> g1 -> Breakpoints -> g2 -> Breakpoints -> CommonPartition g1 g2
mkCommonPartition2 matcher g bg h bh = mkCommonPartition matcher (mkPartition g bg) (mkPartition h bh)

instance (Show g) => Show (Partition g) where
  show pg = combiStr subs_g
    where
      combiStr strs = unwords $ interleavelists strs (replicate (length strs - 1) "|")
      subs_g = map show $ blocks pg

instance (Show g1, Show g2) => Show (CommonPartition g1 g2) where
  show (CGP pg ph) = unlines [show pg, show ph]

writePartition :: CommonPartition GenesIRsR GenesIRsF -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString)
writePartition (CGP pg ph) = (genes_bs_g, irs_bs_g, genes_bs_h, irs_bs_h)
  where
    genes_bs_g = combiBS $ map fst bssg
    irs_bs_g = combiBS $ interleavelists (map snd bssg) ir_breaks_g
    genes_bs_h = combiBS $ map fst bssh
    irs_bs_h = combiBS $ interleavelists (map snd bssh) ir_breaks_h
    combiBS bss = BS.unwords $ interleavelists bss (replicate (length bss - 1) "|")
    ir_breaks_g = map writeIR . breakpointsIR $ pg
    ir_breaks_h = map writeIR . breakpointsIR $ ph
    bssg = map (writeRGenome False) . blocks $ pg
    bssh = map (writeFGenome False) . blocks $ ph

-- Possible position in S of each block from P, indices starting in 0
getBlocksCorrespondence :: (Matcher m g1 g2) => m g1 g2 -> CommonPartition g1 g2 -> [[Int]]
getBlocksCorrespondence matcher (CGP pg ph) = getBlocksCorrespondence_ matcher pg ph

-- Possible position in S of each block from P, indices starting in 0
getBlocksCorrespondence_ :: (Matcher m g1 g2) => m g1 g2 -> Partition g1 -> Partition g2 -> [[Int]]
getBlocksCorrespondence_ matcher pg ph =
  do
    sub_g <- blocks pg
    let sub_hs = blocks ph
    return . map fst . filter (isMatch matcher sub_g . snd) . zip [0 ..] $ sub_hs

data BlockDel = BlockDel
  { b_beg :: Idx,
    b_end :: Idx,
    b_sing :: Bool
  }
  deriving (Show)

instance Eq BlockDel where
  b1 == b2 = (b_beg b1 == b_beg b2) && (b_end b1 == b_end b2)

data CombineWhat = CombineSinSin | CombineSinRep | CombineRepRep

bpsToBlockDels :: (Genome g) => g -> Breakpoints -> Seq BlockDel
bpsToBlockDels g bps = Seq.fromList . zipWith (\a b -> BlockDel (incIdx a) b (any testSingleton [incIdx a .. b])) l_bps $ tail l_bps
  where
    l_bps = 0 : EnumSet.toAscList bps ++ [mkIdx (size g)]
    posMap = positionMap g
    testSingleton i =
      case geneMapLookup (getGene i g) posMap of
        Nothing -> False
        Just pos -> length pos == 1

blockDelsToBps :: Seq BlockDel -> Breakpoints
blockDelsToBps seq_bs = EnumSet.fromList . init $ map (\(BlockDel _ b _) -> b) bs
  where
    bs = toList seq_bs

combine :: (Matcher m g1 g2) => m g1 g2 -> CommonPartition g1 g2 -> CommonPartition g1 g2
combine matcher (CGP (GP g bg) (GP h bh)) = mkCommonPartition2 matcher g bg' h bh'
  where
    (blo_g', blo_h') = combine_ CombineSinSin (blo_g, blo_h)
    blo_g = bpsToBlockDels g bg
    blo_h = bpsToBlockDels h bh
    bg' = blockDelsToBps blo_g'
    bh' = blockDelsToBps blo_h'

    combine_ :: CombineWhat -> (Seq BlockDel, Seq BlockDel) -> (Seq BlockDel, Seq BlockDel)
    combine_ what (bs1, bs2) =
      let (bs1', bs2') = go1 (Empty, bs1, Empty, bs2)
       in if bs1' /= bs1
            then combine_ what (bs1', bs2')
            else case what of
              CombineSinSin -> combine_ CombineSinRep (bs1, bs2)
              CombineSinRep -> combine_ CombineRepRep (bs1, bs2)
              CombineRepRep -> (bs1, bs2)
      where
        -- We have to iterate through the elements of the two sequences,
        -- during the iteration they are moved from the sequence on the right
        -- to the one on the left.
        -- If the sequences are empty there is nothing to do
        go1 (Empty, Empty, Empty, Empty) = (Empty, Empty)
        -- We need at least one more element in the first and second sequence to combine
        go1 (pre, ba1 :<| Empty, Empty, blocks2) = (pre :|> ba1, blocks2)
        go1 (pre1, ba1 :<| bb1 :<| suf1, pre2, ba2 :<| Empty) =
          go1 (pre1 :|> ba1, bb1 :<| suf1, Empty, pre2 :|> ba2)
        -- Now we try to combine
        go1 (pre1, blocks1@(ba1 :<| bb1 :<| suf1), pre2, blocks2@(ba2 :<| bb2 :<| suf2))
          -- test if the two blocks fall into the correct case. If not go to next block in sequence 1.
          | not testReps = go1 (pre1 :|> ba1, bb1 :<| suf1, Empty, pre2 >< blocks2)
          -- test if the block can be combined
          | isMatch matcher g_b1' h_b2' = go1 (pre1, b1' :<| suf1, Empty, (pre2 :|> b2') >< suf2)
          -- next block of sequence 2
          | otherwise = go1 (pre1, blocks1, pre2 :|> ba2, bb2 :<| suf2)
          where
            testReps = case what of
              CombineSinSin -> b_sing ba1 && b_sing bb1
              CombineSinRep -> b_sing ba1 || b_sing bb1
              CombineRepRep -> True

            b1' = BlockDel (b_beg ba1) (b_end bb1) (b_sing ba1 || b_sing bb1)
            b2' = BlockDel (b_beg ba2) (b_end bb2) (b_sing ba2 || b_sing bb2)
            g_b1' = subGenome (b_beg b1') (b_end b1') g
            h_b2' = subGenome (b_beg b2') (b_end b2') h
        go1 (_, _, _, _) = error patternError