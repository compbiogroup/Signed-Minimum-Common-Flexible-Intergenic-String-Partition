{-# LANGUAGE MultiWayIf #-}

module PApprox
  ( approxPartition,
  )
where

import Data.Foldable (foldrM)
import Data.IntMap (IntMap)
import Data.IntMap qualified as IntMap
import Data.List (sortBy)
import Data.Map (Map)
import Data.Map qualified as Map
import Data.Maybe (isJust)
import Data.Set qualified as Set
import Genomes (Genome (size, subGenome), Idx, Matcher (isMatch), idxToInt, mkIdx)
import LocalBase
import Partition (CommonPartition, Partition, blocks, mkCommonPartition, mkPartition, trivialPartition, underlineGenome)

type Correspondents = [Vertex]

data Vertex = Vtx
  { v_size :: Int,
    v_inG :: Bool,
    v_idx :: Idx
  }
  deriving (Eq, Ord)

data BlockMatchGraph g1 g2 = BG (Partition g1) (Partition g2) (Map Vertex (Maybe Int, Correspondents))

data Component = Component
  { c_vtx :: Vertex,
    c_cont_G :: Int,
    c_cont_H :: Int,
    c_perfect :: Bool
  }

approxPartition :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> g1 -> g2 -> CommonPartition g1 g2
approxPartition matcher g h = mkCommonPartition pg ph
  where
    pg_0 = trivialPartition g
    ph_0 = trivialPartition h
    (comps, bmg) = getConnectedComponents (mkBlockMatchGraph matcher pg_0 ph_0)

    (_, BG pg ph _) = until notPart (addBreaks matcher) (comps, bmg)
    notPart = undefined

addBreaks :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> (IntMap Component, BlockMatchGraph g1 g2) -> (IntMap Component, BlockMatchGraph g1 g2)
addBreaks matcher (comps, BG pg ph corres) = (comps', BG pg' ph' corres')
  where
    -- The substrings that may recieve a breakpoint are smallest in components that do not admit a perfect match (this ensures that they are in Tmin)
    cand = takeWhile (\v -> v_size (fst . head $ notPerfs) == v_size (fst v)) notPerfs
    notPerfs = filter isNotPerfect . Map.toAscList $ corres
    isNotPerfect (_, (Nothing, _)) = error patternError
    isNotPerfect (_, (Just c, _)) = not . c_perfect $ comps IntMap.! c
    v = fst . head . sortWith sortKey $ cand
    
    -- We give priority to cut a substring in the genome with more vertices of the
    -- component and break ties with the degree of the vertex.
    sortKey (v_, (Just c_id, neigs)) = (posWeight, dg)
      where
        c = comps IntMap.! c_id
        posWeight =
          if
              | c_cont_G c == c_cont_H c -> 1
              | (c_cont_G c > c_cont_H c) == v_inG v_ -> 0
              | otherwise -> 2
        dg = length neigs
    sortKey (_, (_, _)) = error patternError

    g = underlineGenome pg
    h = underlineGenome ph
    pg' = mkPartition g bg'
    ph' = mkPartition h bh'
    bg' = undefined
    bh' = undefined
    corres' = undefined
    comps' = undefined

-- | Get a set with all connected componets of the Block Match Graph and
-- record the component indices in the Graph
getConnectedComponents :: BlockMatchGraph g1 g2 -> (IntMap Component, BlockMatchGraph g1 g2)
getConnectedComponents (BG pg ph corres) = (comps, BG pg ph corres')
  where
    vertices = Map.assocs corres
    (comps, corres', _) = foldr checkComponent (IntMap.empty, corres, 0) vertices

    checkComponent (v, (maybe_comp, _)) (comps', corres'', comp_id) =
      case maybe_comp of
        Just _ -> (comps', corres'', comp_id)
        Nothing ->
          let (corres''', vs, cont_G, cont_H) = go comp_id v (corres'', [], 0, 0)
              perfect = isJust (testPerfectMatch corres''' vs)
              comp = Component v cont_G cont_H perfect
              comps'' = IntMap.insert comp_id comp comps'
           in (comps'', corres''', comp_id + 1)

    go comp_id v (aux_corres, vs, cont_G, cont_H) =
      let (maybe_comp, neigs) = aux_corres Map.! v
       in case maybe_comp of
            Just _ -> (aux_corres, vs, cont_G, cont_H)
            Nothing ->
              let aux_corres' = Map.insert v (Just comp_id, neigs) aux_corres
                  (cont_G', cont_H') =
                    if v_inG v
                      then (cont_G + 1, cont_H)
                      else (cont_G, cont_H + 1)
               in foldr (go comp_id) (aux_corres', v : vs, cont_G', cont_H') neigs

    testPerfectMatch aux_corres = foldrM (selectPBlock Set.empty) Map.empty
      where
        selectPBlock seen v blocksFromG =
          if perfect then Just blocksFromG' else Nothing
          where
            (_, blocksFromG', perfect) = selectPBlock' seen v blocksFromG

        selectPBlock' seen v = testCandidates v neigs seen
          where
            (_, neigs) = aux_corres Map.! v

        testCandidates _ [] seen blocksFromG = (seen, blocksFromG, False)
        testCandidates v (u : us) seen blocksFromG =
          if
              | u `Set.member` seen -> testCandidates v us seen blocksFromG
              | u `Map.notMember` blocksFromG -> (seen', blocksFromG', True)
              | perfect -> (seen_rec, blocksFromG_rec', True)
              | otherwise -> testCandidates v us seen_rec blocksFromG_rec
          where
            seen' = Set.insert u seen
            blocksFromG' = Map.insert u v blocksFromG
            (seen_rec, blocksFromG_rec, perfect) = selectPBlock' seen' (blocksFromG Map.! u) blocksFromG
            blocksFromG_rec' = Map.insert u v blocksFromG_rec

mkBlockMatchGraph :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> Partition g1 -> Partition g2 -> BlockMatchGraph g1 g2
mkBlockMatchGraph matcher pg ph = BG pg ph corres
  where
    bsG = blocks pg
    bsH = blocks ph
    acc_block_sizes :: (Genome g) => [g] -> [(g, Int)]
    acc_block_sizes bsK = zip bsK . scanl (+) 0 . map size $ bsK
    corres = Map.fromList $ corresList bsG getCorresInH True ++ corresList bsH getCorresInG False
    getCorresInH = getCorresIn bsH (isMatch matcher) False
    getCorresInG = getCorresIn bsG (flip (isMatch matcher)) True

    corresList bsQ getCorresInK inG = do
      (b, prev_sizes) <- acc_block_sizes bsQ
      beg <- [1 .. mkIdx (size b)]
      end <- [beg + 1 .. mkIdx (size b)]
      let sub = subGenome beg end b
      return (Vtx (size sub) inG (mkIdx prev_sizes + beg), (Nothing, getCorresInK sub))

    getCorresIn bsK test inG sub = do
      (b, prev_sizes) <- acc_block_sizes bsK
      beg <- [1 .. mkIdx (size b)]
      let end = beg + mkIdx (size sub)
      [ Vtx (size sub) inG (mkIdx prev_sizes + beg)
        | end < mkIdx (size b) && test sub (subGenome beg end b)
        ]