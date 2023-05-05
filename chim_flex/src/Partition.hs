{-# LANGUAGE ExplicitForAll #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}

module Partition
  ( greedyPart,
    soarPartition,
    getPartition,
    writePartition,
    getBlocksMatchGraph,
    longestSubstring,
    commonPrefix,
    suboptimalRuleInterval,
    suboptimalRulePairs,
  )
where

import Control.Applicative ((<|>))
import Data.Bits (Bits, complement, setBit, testBit, zeroBits)
import Data.ByteString.Char8 qualified as BS
import Data.Foldable (toList)
import Data.IntSet (IntSet)
import Data.IntSet qualified as IntSet
import Data.List qualified as List
import Data.Maybe (catMaybes, fromMaybe, listToMaybe, mapMaybe, maybeToList)
import Genomes (Gene, GeneMap, GenesIRsF, GenesIRsR, Genome (..), Idx, IntergenicGenome (..), Matcher (..), RigidFlexibleReverseMatcher (..), decIdx, geneMapLookup, idxDist, idxToInt, incIdx, mkIdx, positionMap, writeFGenome, writeIR, writeRGenome)
import LocalBase
import Text.Printf (PrintfArg, printf)

type Breakpoint = [Idx]

data Partition g1 g2 where
  GenomePartition :: (Genome g1, Genome g2) => g1 -> g2 -> Breakpoint -> Breakpoint -> Partition g1 g2

instance Show (Partition GenesIRsR GenesIRsF) where
  show (GenomePartition g h bg bh) = unlines [combiStr subs_g, combiStr subs_h]
    where
      combiStr strs = unwords $ interleavelists strs (replicate (length strs - 1) "|")
      subs_g = zipWith (getSub show g) bg $ tail bg
      subs_h = zipWith (getSub show h) bh $ tail bh
      getSub write g' i succi = write $ subGenome (incIdx i) succi g'

getPartition :: GenesIRsR -> GenesIRsF -> Partition GenesIRsR GenesIRsF
getPartition = soarPartition RFRM

writePartition :: Partition GenesIRsR GenesIRsF -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString)
writePartition (GenomePartition g h bg bh) = (genes_bs_g, irs_bs_g, genes_bs_h, irs_bs_h)
  where
    genes_bs_g = combiBS $ map fst bssg
    irs_bs_g = combiBS $ interleavelists (map snd bssg) ir_breaks_g
    genes_bs_h = combiBS $ map fst bssh
    irs_bs_h = combiBS $ interleavelists (map snd bssh) ir_breaks_h
    combiBS bss = BS.unwords $ interleavelists bss (replicate (length bss - 1) "|")
    ir_breaks_g = map (writeIR . (`getIR` g)) (init . tail $ bg)
    ir_breaks_h = map (writeIR . (`getIR` h)) (init . tail $ bh)
    bssg = zipWith (getSub (writeRGenome False) g) bg $ tail bg
    bssh = zipWith (getSub (writeFGenome False) h) bh $ tail bh
    getSub write g' i succi = write $ subGenome (incIdx i) succi g'

getBlocksMatchGraph :: RigidFlexibleReverseMatcher GenesIRsR GenesIRsF -> Partition GenesIRsR GenesIRsF -> [[Int]]
getBlocksMatchGraph macher (GenomePartition g h bg bh) =
  do
    sub_g <- zipWith (getSub g) bg $ tail bg
    let sub_hs = zipWith (getSub h) bh $ tail bh
    return . map fst . filter (isMatch macher sub_g . snd) . zip [0 ..] $ sub_hs
  where
    getSub g' i succi = subGenome (incIdx i) succi g'

-- | Greedy algorithm for string partition with intergenic regions or
-- flexible intergeic regions. If withSingleton is True it looks first
-- for matched blocks that contain a singleton.
greedyPart :: (Genome g1, Genome g2, Matcher m g1 g2) => Bool -> m g1 g2 -> g1 -> g2 -> Partition g1 g2
greedyPart withSingleton matcher g h = GenomePartition g h (cleanList final_breaksG) (cleanList final_breaksH)
  where
    cleanList = map head . List.group . List.sort
    (final_breaksG, final_breaksH) = getAllLS all_singletons IntSet.empty IntSet.empty [] []
    getAllLS singletons genesOnBlocksG genesOnBlocksH breaksG breaksH =
      case longestSubstring' singletons of
        (Nothing, Nothing) -> (breaksG, breaksH)
        (Just other_singletons, Nothing) -> getAllLS other_singletons genesOnBlocksG genesOnBlocksH breaksG breaksH
        (maybe_other_singletons, Just (((g_beg, g_end), (h_beg, h_end)), genesOnBlocksG', genesOnBlocksH')) ->
          let breaksG' = (decIdx g_beg : g_end : breaksG)
              breaksH' = (decIdx h_beg : h_end : breaksH)
           in getAllLS (fromMaybe [] maybe_other_singletons) genesOnBlocksG' genesOnBlocksH' breaksG' breaksH'
      where
        longestSubstring' (singleton : other_singletons) = (Just other_singletons, longestSubstring (Just singleton) matcher genesOnBlocksG genesOnBlocksH g h)
        longestSubstring' [] = (Nothing, longestSubstring Nothing matcher genesOnBlocksG genesOnBlocksH g h)

    all_singletons =
      if withSingleton
        then filter (singletonOnBoth posMapG posMapH) . toList $ alphabet g
        else []
    posMapG = positionMap g
    posMapH = positionMap h

singletonOnBoth :: GeneMap [Idx] -> GeneMap [Idx] -> Gene -> Bool
singletonOnBoth posMapG posMapH gene =
  (case geneMapLookup gene posMapG of Nothing -> False; Just pos -> length pos == 1)
    && (case geneMapLookup gene posMapH of Nothing -> False; Just pos -> length pos == 1)

longestSubstring :: forall m g1 g2. (Genome g1, Genome g2, Matcher m g1 g2) => Maybe Gene -> m g1 g2 -> IntSet -> IntSet -> g1 -> g2 -> Maybe (((Idx, Idx), (Idx, Idx)), IntSet, IntSet)
-- Find longest substring of l1 and l2, ignore genes already in a block
-- returns the indices of the correspondent subgenomes in g and h (first and last gene).
-- If singleton has a value (Just x), only consider blocks containing that gene x.
longestSubstring maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h =
  case maybe_inds of
    Nothing -> Nothing
    Just inds@((g_beg, g_end), (h_beg, h_end)) ->
      Just
        ( inds,
          foldr (IntSet.insert . idxToInt) genesOnBlocksG [g_beg .. g_end],
          foldr (IntSet.insert . idxToInt) genesOnBlocksH [h_beg .. h_end]
        )
  where
    maybe_inds =
      maxWith (uncurry idxDist . fst)
        . mapMaybe (longestMatch . mkIdx)
        . filter (`IntSet.notMember` genesOnBlocksG)
        $ [1 .. size g]

    longestMatch :: Idx -> Maybe ((Idx, Idx), (Idx, Idx))
    -- find longest match between a subgenome of h and a suffix of g,
    -- ignore genes already in a block,
    -- returns the indices of the correspondent subgenomes in g and h (first and last gene)
    longestMatch ig = maxWith (uncurry idxDist . fst) $ do
      ih <- map mkIdx . filter (`IntSet.notMember` genesOnBlocksH) $ [1 .. size h]
      maybeToList $ commonPrefix maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h ig ih

commonPrefix :: (Matcher m g1 g2, Genome g1, Genome g2) => Maybe Gene -> m g1 g2 -> IntSet -> IntSet -> g1 -> g2 -> Idx -> Idx -> Maybe ((Idx, Idx), (Idx, Idx))
-- find common prefix of subgenome of g starting with ig
-- and subgenome of h starting with ih
-- returns Nothing if the prefix is empty
-- returns indices for the ends of both subgenomes otherwise
commonPrefix maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h ig ih = do
  (endCPG, endCPH) <-
    if
        | isDirectMatch matcher (subGenome ig ig g) (subGenome ih ih h) -> Just $ commonPrefix' False ig ih
        | isReverseMatch matcher (subGenome ig ig g) (subGenome ih ih h) -> Just $ commonPrefix' True ig ih
        | otherwise -> Nothing
  case maybe_singleton of
    Just singleton ->
      if singleton `isGene` subGenome ig endCPG g
        then return ((ig, endCPG), (ih, endCPH))
        else Nothing
    Nothing -> return ((ig, endCPG), (ih, endCPH))
  where
    commonPrefix' :: Bool -> Idx -> Idx -> (Idx, Idx)
    commonPrefix' rev jg jh =
      if
          | jg == mkIdx (size g) -> (jg, jh)
          | jh == mkIdx (size h) -> (jg, jh)
          | idxToInt jg' `IntSet.member` genesOnBlocksG -> (jg, jh)
          | idxToInt jh' `IntSet.member` genesOnBlocksH -> (jg, jh)
          | not (testMatch matcher (subGenome jg jg' g) (subGenome jh jh' h)) -> (jg, jh)
          | otherwise -> commonPrefix' rev jg' jh'
      where
        jg' = incIdx jg
        jh' = incIdx jh
        testMatch = if rev then isReverseMatch else isDirectMatch

combine :: Partition g1 g2 -> Partition g1 g2
combine (GenomePartition g h bg bh) = undefined

soarPartition :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> g1 -> g2 -> Partition g1 g2
soarPartition matcher g_ h_ = part -- combine part
  where
    (g, h) = uncurry (suboptimalRulePairs matcher) (suboptimalRuleInterval matcher g_ h_)
    nG = size g
    nH = size h
    part = GenomePartition g h (cleanList breaksG) (cleanList breaksH)
    cleanList = map head . List.group . List.sort
    graph = makePMGraph4 matcher g h
    (breaksG, breaksH) = traceValue . (getBpsFromIS nG nH graph) . independentSet $ graph

newtype BitMask = BM Integer deriving (Eq, Ord, Bits, PrintfArg)

instance Show BitMask where
  show = printf "%llb"

class PMGraph pmg where
  independentSet :: pmg -> BitMask
  vertexCover :: pmg -> BitMask

  independentSet = complement . vertexCover
  vertexCover = complement . independentSet

type Vertex = Int

type Edge = (Vertex, Vertex)

type PairMatch g1 g2 = (g1, g2, Idx, Idx, Vertex)

data PMGraph4 g1 g2 = PMGraph4 [Edge] [PairMatch g1 g2]

instance (Show g1, Show g2) => Show (PMGraph4 g1 g2) where
  show (PMGraph4 edges pairMatches) = show edges ++ "\n" ++ show pairMatches

isConflict :: (Matcher m g1 g2) => m g1 g2 -> PairMatch g1 g2 -> PairMatch g1 g2 -> Bool
isConflict matcher pm@(due1, due2, idx1, idx2, _) pm'@(_, _, idx1', idx2', _) =
  if idx1 > idx1'
    then isConflict matcher pm' pm
    else
      idx1 == idx1' && (idx2 /= idx2')
        || idx1 + 1 < idx1' && (List.intersect [idx2, idx2 + 1] [idx2', idx2' + 1] /= [])
        || idx1 + 1 == idx1'
          && not
            ( idx2 + 1 == idx2' && isDirectMatch matcher due1 due2
                || idx2 == idx2' + 1 && isReverseMatch matcher due1 due2
            )

makePMGraph4 :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> g1 -> g2 -> PMGraph4 g1 g2
makePMGraph4 matcher g h = PMGraph4 edges pms
  where
    edges = catMaybes (testEdge <$> pms <*> pms)
    testEdge pm1@(_, _, _, _, v) pm2@(_, _, _, _, u) =
      if isConflict matcher pm1 pm2
        then Just (v, u)
        else Nothing
    pms = zipWith addVertex [1 ..] . catMaybes $ (check <$> duesG <*> duesH)
    addVertex v (due1, due2, idx1, idx2) = (due1, due2, mkIdx idx1, mkIdx idx2, v)
    check (idx1, due1) (idx2, due2) =
      if isMatch matcher due1 due2
        then Just (due1, due2, idx1, idx2)
        else Nothing
    duesG = map (\i -> (i, subGenome (mkIdx i) (mkIdx (i + 1)) g)) [1 .. (size g - 1)]
    duesH = map (\i -> (i, subGenome (mkIdx i) (mkIdx (i + 1)) h)) [1 .. (size h - 1)]

instance PMGraph (PMGraph4 g1 g2) where
  vertexCover (PMGraph4 edges _) = foldr coverEdges zeroBits edges
    where
      coverEdges :: Edge -> BitMask -> BitMask
      coverEdges (u, v) onCover
        | not (testBit onCover u) && not (testBit onCover v) = setBit (setBit onCover u) v
        | otherwise = onCover

getBpsFromIS :: Int -> Int -> PMGraph4 g1 g2 -> BitMask -> ([Idx], [Idx])
getBpsFromIS nG nH (PMGraph4 _ pms) indSet = (bpsG, bpsH)
  where
    bpsG = [1 .. mkIdx nG] List.\\ getDues False
    bpsH = [1 .. mkIdx nH] List.\\ getDues True
    getDues isSecond =
      mapMaybe
        ( \(_, _, idx1, idx2, v) ->
            if testBit indSet v
              then Just (if isSecond then idx2 else idx1)
              else Nothing
        )
        pms

-- | Match two correspondent intervals a^l...a^r from g and b^l...b^r from h.
suboptimalRuleInterval :: (Matcher m g1 g2, Genome g1, Genome g2) => m g1 g2 -> g1 -> g2 -> (g1, g2)
suboptimalRuleInterval matcher g h = foldr checkAndReplace (g, h) is_js
  where
    posMapG = positionMap g
    posMapH = positionMap h
    originallySingleton i = singletonOnBoth posMapG posMapH (getGene i g)

    -- Check possible interval matches and replace than with singletons
    checkAndReplace (i, j) (g', h') =
      case checkDirMatch i j <|> checkRevMatch i j <|> checkMiddleMatch i j of
        Nothing -> (g', h')
        Just (al, ar, bl, br, isRev) ->
          let (g'', singletons_g) = makeSingletons h [al + 1 .. ar - 1] g'
              singletons_h =
                if isRev
                  then reverse . map (invGene g) $ singletons_g
                  else singletons_g
              h'' = setGenes [bl + 1 .. br - 1] singletons_h h'
           in (g'', h'')

    -- candidates for indices of a^l in a^l...a^r
    is = filter originallySingleton [1 .. mkIdx (size g)]
    -- combine is with the correspondent indices in h
    is_js =
      map
        ( \i -> case geneMapLookup (getGene i g) posMapH of
            Just [j] -> (i, j)
            _ -> error patternError
        )
        is

    -- check if the interval is a direct match (a^l...a^r = b^l...b^r) and
    -- returns indices of (a^l,a^r,b^l,b^r) and the value False indicating that
    -- the match is not reverse.
    checkDirMatch i j = do
      k <- k_result
      if k > 1
        then Just (i, i + k, j, j + k, False)
        else Nothing
      where
        -- k is the displacement from a^l to a^r
        ks = [1 .. (min (mkIdx (size g) - i) (mkIdx (size h) - j))]
        matches = map fst . takeWhile snd $ map (\k -> (k, isDirectMatch matcher (subGenome (i + k - 1) (i + k) g) (subGenome (j + k - 1) (j + k) h))) ks
        k_result =
          (\k -> if k > 1 then Just k else Nothing)
            =<< List.find (\k -> originallySingleton (i + k)) matches

    -- check if the interval is a reverse match (a^l...a^r = -b^r...-b^l) and
    -- returns indices of (a^l,a^r,b^l,b^r) and the value True indicating that
    -- the match is reverse.
    checkRevMatch i j = do
      k <- k_result
      if k > 1
        then Just (i, i + k, j - k, j, True)
        else Nothing
      where
        -- k is the displacement from a^l to a^r
        ks = [1 .. (min (mkIdx (size g) - i) (j - 1))]
        matches = map fst . takeWhile snd $ map (\k -> (k, isReverseMatch matcher (subGenome (i + k - 1) (i + k) g) (subGenome (j - k) (j - k + 1) h))) ks
        k_result =
          (\k -> if k > 1 then Just k else Nothing)
            =<< List.find (\k -> originallySingleton (i + k)) matches

    -- check if the middle is a reverse match
    -- (a^la^{l+1}...a^{r-1}a^r = b^l-b^{r-1}...-b^{l+1}b^r) and
    -- returns indices of (a^l,a^r,b^l,b^r) and the value True indicating that
    -- the match is reverse.
    checkMiddleMatch i j = do
      idxAr <- getIdxAr
      k <- k_result idxAr
      Just (i, idxAr, j, j + k, True)
      where
        getIdxAr = List.find originallySingleton [i + 1 .. mkIdx (size g)]
        -- k is the displacement from a^l to a^r
        ks idxAr = [2 .. idxAr - i - 1]
        matches idxAr = map fst . takeWhile snd . map (\k -> (k, isReverseMatch matcher (subGenome (idxAr - k) (idxAr - k + 1) g) (subGenome (j + k - 1) (j + k) h))) $ ks idxAr
        k_result idxAr =
          if idxAr - i > mkIdx (size h) - j
            || not
              ( isDirectMatch
                  matcher
                  (subGenome i i g)
                  (subGenome j j h)
              )
            || not
              ( isDirectMatch
                  matcher
                  (subGenome idxAr idxAr g)
                  (subGenome (j + idxAr - i) (j + idxAr - i) h)
              )
            then Nothing
            else listToMaybe (matches idxAr)

-- | Match pair of a singleton and a replica (a^ha^t).

--- The replica in mapped into a singleton.
suboptimalRulePairs :: (Matcher m g1 g2, Genome g1, Genome g2) => m g1 g2 -> g1 -> g2 -> (g1, g2)
suboptimalRulePairs matcher g h = foldr replaceWithSingletons (g, h) (match_ah ++ match_at)
  where
    posMapG = positionMap g
    posMapH = positionMap h
    originallySingleton i = singletonOnBoth posMapG posMapH (getGene i g)

    -- replace replica of each matched pair with a singleton
    replaceWithSingletons (a, b, isRev) (g', h') =
      let (g'', singletons_g) = makeSingletons h [a] g'
          singletons_h =
            if isRev
              then map (invGene g) singletons_g
              else singletons_g
          h'' = setGenes [b] singletons_h h'
       in (g'', h'')

    -- candidates for indices of singletons
    is = filter originallySingleton [1 .. mkIdx (size g)]

    -- a^h is a singleton in both genomes and a^t is a replica
    is_ah = filter (\i -> i < mkIdx (size g) && testReplicaG (i + 1)) is

    -- find indices from G and H that must be turn into singletons in this case
    match_ah = mapMaybe (uncurry checkPair_ah) $ findInH is_ah

    checkPair_ah i j =
      if
          | j < mkIdx (size h) && isDirectMatch matcher (subGenome i (i + 1) g) (subGenome j (j + 1) h) -> Just (i + 1, j + 1, False)
          | j > 1 && isReverseMatch matcher (subGenome i (i + 1) g) (subGenome (j - 1) j h) -> Just (i + 1, j - 1, True)
          | otherwise -> Nothing

    -- a^t is a singleton in both genomes and a^h is a replica
    is_at = filter (\i -> i > 1 && testReplicaG (i - 1)) is

    -- find indices from G and H that must be turn into singletons in this case
    match_at = mapMaybe (uncurry checkPair_at) $ findInH is_at

    checkPair_at i j =
      if
          | j < mkIdx (size h) && isReverseMatch matcher (subGenome (i - 1) i g) (subGenome j (j + 1) h) -> Just (i - 1, j + 1, True)
          | j > 1 && isDirectMatch matcher (subGenome (i - 1) i g) (subGenome (j - 1) j h) -> Just (i - 1, j - 1, False)
          | otherwise -> Nothing

    -- combine is with the correspondent indices in h
    findInH =
      map
        ( \i -> case geneMapLookup (getGene i g) posMapH of
            Just [j] -> (i, j)
            _ -> error patternError
        )

    -- test if gene at position i is a replica in the genome g
    testReplicaG i =
      case geneMapLookup (getGene i g) posMapG of
        Nothing -> False
        Just pos -> length pos > 1
