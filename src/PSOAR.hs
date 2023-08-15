{-# LANGUAGE GeneralisedNewtypeDeriving #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}

-- \|
-- Module      : PGreedy
-- Description : Implementation of the algorithm in the SOAR framework for string
-- partition adapted to include intergenic region information.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module PSOAR
  ( soarPartition,
    suboptimalRuleInterval,
    suboptimalRulePairs,
  )
where

import Control.Applicative ((<|>))
import Data.Bits (Bits (complement, setBit, testBit, zeroBits))
import Data.EnumSet (EnumSet)
import Data.EnumSet qualified as EnumSet
import Data.List qualified as List
import Data.Maybe (catMaybes, listToMaybe, mapMaybe)
import Genomes (Genome (..), Idx, Matcher (..), geneMapLookup, mkIdx, positionMap, singletonOnBoth)
import LocalBase
import Partition (CommonPartition, combine, mkCommonPartition2)
import Text.Printf (PrintfArg, printf)

soarPartition :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> g1 -> g2 -> CommonPartition g1 g2
soarPartition matcher g_ h_ = part -- combine matcher part
  where
    (g, h) = uncurry (suboptimalRulePairs matcher) (suboptimalRuleInterval matcher g_ h_)
    -- replace the previous line with this next one to disable the use of
    -- the suboptimal rules
    -- (g, h) = (g_, h_)
    nG = size g
    nH = size h
    part = mkCommonPartition2 matcher g breaksG h breaksH
    graph = makePMGraph4 matcher g h
    (breaksG, breaksH) = getBpsFromIS nG nH graph . independentSet $ graph

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

getBpsFromIS :: Int -> Int -> PMGraph4 g1 g2 -> BitMask -> (EnumSet Idx, EnumSet Idx)
getBpsFromIS nG nH (PMGraph4 _ pms) indSet = (bpsG, bpsH)
  where
    bpsG = EnumSet.fromList $ [1 .. mkIdx nG - 1] List.\\ getDues False
    bpsH = EnumSet.fromList $ [1 .. mkIdx nH - 1] List.\\ getDues True
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
      case checkDirMatch i j g' h' <|> checkRevMatch i j g' h' <|> checkMiddleMatch i j g' h' of
        Nothing -> (g', h')
        Just (al, ar, bl, br, isRev) ->
          let (g'', singletons_g) = makeSingletons h' [al + 1 .. ar - 1] g'
              singletons_h =
                if isRev
                  then reverse . map (invGene g') $ singletons_g
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
    checkDirMatch i j g' h' = do
      k <- k_result
      if k > 1
        then Just (i, i + k, j, j + k, False)
        else Nothing
      where
        -- k is the displacement from a^l to a^r
        ks = [1 .. (min (mkIdx (size g') - i) (mkIdx (size h') - j))]
        matches = map fst . takeWhile snd $ map (\k -> (k, isDirectMatch matcher (subGenome (i + k - 1) (i + k) g') (subGenome (j + k - 1) (j + k) h'))) ks
        k_result =
          (\k -> if k > 1 then Just k else Nothing)
            =<< List.find (\k -> originallySingleton (i + k)) matches

    -- check if the interval is a reverse match (a^l...a^r = -b^r...-b^l) and
    -- returns indices of (a^l,a^r,b^l,b^r) and the value True indicating that
    -- the match is reverse.
    checkRevMatch i j g' h' = do
      k <- k_result
      if k > 1
        then Just (i, i + k, j - k, j, True)
        else Nothing
      where
        -- k is the displacement from a^l to a^r
        ks = [1 .. (min (mkIdx (size g') - i) (j - 1))]
        matches = map fst . takeWhile snd $ map (\k -> (k, isReverseMatch matcher (subGenome (i + k - 1) (i + k) g') (subGenome (j - k) (j - k + 1) h'))) ks
        k_result =
          (\k -> if k > 1 then Just k else Nothing)
            =<< List.find (\k -> originallySingleton (i + k)) matches

    -- check if the middle is a reverse match
    -- (a^la^{l+1}...a^{r-1}a^r = b^l-b^{r-1}...-b^{l+1}b^r) and
    -- returns indices of (a^l,a^r,b^l,b^r) and the value True indicating that
    -- the match is reverse.
    checkMiddleMatch i j g' h' = do
      k <- listToMaybe $ mapMaybe k_result getIdxAr
      Just (i, i + k, j, j + k, True)
      where
        getIdxAr = filter originallySingleton [i + 2 .. mkIdx (size g')]
        -- k is the displacement from a^l to a^r
        ks idxAr = [2 .. idxAr - i - 1]
        matches idxAr = all (\k -> isReverseMatch matcher (subGenome (idxAr - k) (idxAr - k + 1) g') (subGenome (j + k - 1) (j + k) h')) $ ks idxAr
        k_result idxAr =
          if idxAr - i > mkIdx (size h') - j
            || not
              ( isDirectMatch
                  matcher
                  (subGenome i i g')
                  (subGenome j j h')
              )
            || not
              ( isDirectMatch
                  matcher
                  (subGenome idxAr idxAr g')
                  (subGenome (j + idxAr - i) (j + idxAr - i) h')
              )
            || idxAr - i == 2
              && not
                ( isReverseMatch
                    matcher
                    (subGenome (i + 1) (i + 1) g')
                    (subGenome (j + 1) (j + 1) h')
                )
            || not (matches idxAr)
            then Nothing
            else Just (idxAr - i)

-- | Match pair of a singleton and a replica (a^ha^t).

--- The replica in mapped into a singleton.
suboptimalRulePairs :: (Matcher m g1 g2, Genome g1, Genome g2) => m g1 g2 -> g1 -> g2 -> (g1, g2)
suboptimalRulePairs matcher g h = foldr checkAndReplace (g, h) (match_ah ++ match_at)
  where
    posMapG = positionMap g
    posMapH = positionMap h
    originallySingleton i = singletonOnBoth posMapG posMapH (getGene i g)

    -- replace replica of each matched pair with a singleton
    checkAndReplace pairsToCheck (g', h') =
      case pairsToCheck g' h' of
        Nothing -> (g', h')
        Just (a, b, isRev) ->
          let (g'', singletons_g) = makeSingletons h' [a] g'
              singletons_h =
                if isRev
                  then map (invGene g') singletons_g
                  else singletons_g
              h'' = setGenes [b] singletons_h h'
           in (g'', h'')

    -- candidates for indices of singletons
    is = filter originallySingleton [1 .. mkIdx (size g)]

    -- a^h is a singleton in both genomes and a^t is a replica
    is_ah = filter (\i -> i < mkIdx (size g) && testReplicaG (i + 1)) is

    -- find indices from G and H that must be turn into singletons in this case
    match_ah = map (\(i, j) g' h' -> checkPair_ah g' h' i j) $ findInH is_ah

    checkPair_ah g' h' i j =
      if
          | j < mkIdx (size h') && isDirectMatch matcher (subGenome i (i + 1) g') (subGenome j (j + 1) h') -> Just (i + 1, j + 1, False)
          | j > 1 && isReverseMatch matcher (subGenome i (i + 1) g') (subGenome (j - 1) j h') -> Just (i + 1, j - 1, True)
          | otherwise -> Nothing

    -- a^t is a singleton in both genomes and a^h is a replica
    is_at = filter (\i -> i > 1 && testReplicaG (i - 1)) is

    -- find indices from G and H that must be turn into singletons in this case
    match_at = map (\(i, j) g' h' -> checkPair_at g' h' i j) $ findInH is_at

    checkPair_at g' h' i j =
      if
          | j < mkIdx (size h') && isReverseMatch matcher (subGenome (i - 1) i g') (subGenome j (j + 1) h') -> Just (i - 1, j + 1, True)
          | j > 1 && isDirectMatch matcher (subGenome (i - 1) i g') (subGenome (j - 1) j h') -> Just (i - 1, j - 1, False)
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
