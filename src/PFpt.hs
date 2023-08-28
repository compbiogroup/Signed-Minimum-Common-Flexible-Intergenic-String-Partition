{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TupleSections #-}

-- \|
-- Module      : Genomes
-- Description : Algorithm for genome partition with an polinomial complexity
-- parameterized by the partition size, maximum occurrence and/or number of distinct
-- elements between the alphabets. Currently the sample graph used in this algorithm
-- does not work for an unsigned Reverse partition. Consequently, the algorithm should
-- be used only for direct and signed partition problems.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module PFpt
  ( fptPartition,
  )
where

import Control.Concurrent (newMVar, readMVar, swapMVar)
import Control.Concurrent.Async (race)
import Control.Concurrent.MVar (MVar)
import Control.Exception (assert)
import Control.Monad (void)
import Data.Either (isRight)
import Data.EnumSet qualified as EnumSet
import Data.Map qualified as Map
import Data.Maybe (isNothing, listToMaybe)
import Data.MultiMap (MultiMap)
import Data.MultiMap qualified as MMap
import Data.Set qualified as Set
import GHC.Conc (threadDelay)
import Genomes (Genome (getGene, size, subGenome), Idx, Matcher (isDirectMatch, isReverseMatch), decIdx, incIdx, mkIdx)
import LocalBase
import Partition (CommonPartition, mkCommonPartition2)

data Vertex = Vtx
  { vInG :: Bool,
    vIdx :: Idx
  }
  deriving (Eq, Ord, Show)

-- In all Idx Pairs in Edge and Samples we assume that the first index is
-- from G and the second from H
data Edge = Edge
  { eOther :: Vertex,
    eVpair :: (Idx, Idx),
    eParent :: BlockWitness,
    eBlack :: Bool
  }
  deriving (Eq, Ord, Show)

type Sample = [BlockWitness]

data BlockWitness = BW
  { w_left_pair :: (Idx, Idx), -- leftmost pair (looking in G) parallel to black_pair
    w_black_pair :: (Idx, Idx),
    w_right_pair :: (Idx, Idx), -- rightmost pair (looking in G) parallel to black_pair
    w_inverted :: Bool
  }
  deriving (Eq, Ord, Show)

data SampleGraph m g1 g2 = SG
  { sg_matcher :: m g1 g2,
    sg_g :: g1,
    sg_h :: g2,
    sample :: Sample,
    edges :: MultiMap Vertex Edge
  }

instance (Show g1, Show g2) => Show (SampleGraph m g1 g2) where
  show SG {..} = unlines [show sg_g, show sg_h, show sample, edges_str]
    where
      edges_str = unlines . map show . MMap.toList $ edges

-- TODO: To deal with unbalance strings the branching rule 2 still have to be
-- implemented and a text must be include on the rare odd path to ensure that they
-- are rare
fptPartition :: (Genome g1, Genome g2, Matcher m g1 g2) => Int -> m g1 g2 -> g1 -> g2 -> IO (Maybe (CommonPartition g1 g2), Bool)
fptPartition time_limit matcher g h = do
  best_sg <- newMVar Nothing
  full_comp <- isRight <$> race (threadDelay time_limit) (runFpt best_sg sg0)
  maybe_sg <- readMVar best_sg
  return
    . (,full_comp)
    . fmap
      ( sampleGraphToPartition
          . (\sg' -> assert (checkDegrees 1 1 sg') sg')
          . finalizeCorrespondence
          . (\sg' -> assert (checkDegrees 1 2 sg') sg')
      )
    $ maybe_sg
  where
    sg0 = (\sg -> assert (checkDegrees 0 2 sg) sg) $ createSampleGraph matcher g h sample0
    sample0 = [] -- TODO: add an initial sample with the singletons of both strings

runFpt :: (Genome g1, Genome g2, Matcher m g1 g2) => MVar (Maybe (SampleGraph m g1 g2)) -> SampleGraph m g1 g2 -> IO ()
runFpt best_sg = go
  where
    go !sg = do
      test <- toBig sg
      if
          | test -> return ()
          | isNothing maybe_oddPath -> selectSG sg
          | null branch_edges -> return ()
          | otherwise ->
              mapM_
                ( \edge ->
                    let sg' = addBlackEdge edge sg
                     in assert (checkDegrees 0 2 sg') $ go sg'
                )
                branch_edges
      where
        branch_edges = maybe [] (edgesFromRareOddPath sg) maybe_oddPath
        maybe_oddPath = getRareOddPath sg

    toBig sg = do
      maybe_best_sg <- readMVar best_sg
      case maybe_best_sg of
        Nothing -> return False
        Just best_sg' -> return $ length (sample sg) >= length (sample best_sg')
    selectSG sg = do
      maybe_best_sg <- readMVar best_sg
      case maybe_best_sg of
        Nothing ->
          void (swapMVar best_sg (Just sg))
        Just best_sg' -> do
          _ <- swapMVar best_sg (if length (sample sg) < length (sample best_sg') then Just sg else Just best_sg')
          return ()

sampleGraphToPartition :: (Matcher m g1 g2, Genome g1, Genome g2) => SampleGraph m g1 g2 -> CommonPartition g1 g2
sampleGraphToPartition (SG {..}) = mkCommonPartition2 sg_matcher sg_g bg sg_h bh
  where
    bg = foldr (processBlackEdges True) EnumSet.empty [2 .. mkIdx (size sg_g)]
    bh = foldr (processBlackEdges False) EnumSet.empty [2 .. mkIdx (size sg_h)]
    processBlackEdges inG idx breaks =
      if eParent e1 /= eParent e2
        then EnumSet.insert (idx - 1) breaks
        else breaks
      where
        e1 = head . (\l -> if null l then undefined else l) $ edges MMap.! Vtx inG (idx - 1)
        e2 = head . (\l -> if null l then undefined else l) $ edges MMap.! Vtx inG idx

-- return vertices from a rare odd path if such path exists
getRareOddPath :: (Genome g1) => SampleGraph m g1 g2 -> Maybe [Vertex]
getRareOddPath = listToMaybe . getRareOddPaths

-- return vertices from each rare odd path
getRareOddPaths :: (Genome g1) => SampleGraph m g1 g2 -> [[Vertex]]
getRareOddPaths (SG {..}) = result_list
  where
    result_list =
      snd $
        foldr (getOddPathWithVertex . Vtx True) (Set.empty, []) [1 .. mkIdx (size sg_g)]
    getOddPathWithVertex v (seen, vs) =
      if
          | v `Set.member` seen -> (seen, vs)
          | isOdd -> (seen', vertices : vs)
          | otherwise -> (seen', vs)
      where
        isOdd = length vertices `mod` 2 == 1
        (vertices, seen') =
          case edges MMap.! v of
            [] -> ([v], aux_seen)
            [e] ->
              let (verticesRest, aux_seen') = getPath (eOther e) aux_seen
               in (v : verticesRest, aux_seen')
            [e1, e2] ->
              let (vertices1, aux_seen') = getPath (eOther e1) aux_seen
                  (vertices2, aux_seen'') = getPath (eOther e2) aux_seen'
               in (vertices1 ++ (v : vertices2), aux_seen'')
            _ -> error patternError
          where
            aux_seen = Set.insert v seen

    getPath v aux_seen =
      if v `Set.member` aux_seen
        then ([], aux_seen)
        else case filter ((`Set.notMember` aux_seen) . eOther) (edges MMap.! v) of
          [] -> ([v], aux_seen')
          [e] ->
            let (verticesRest, aux_seen'') = getPath (eOther e) aux_seen'
             in (v : verticesRest, aux_seen'')
          _ -> error patternError
      where
        aux_seen' = Set.insert v aux_seen

-- return true if the black edges from the witnesses share and edge
checkConflict :: BlockWitness -> BlockWitness -> Bool
checkConflict bw bw' = idxG == idxG' || idxH == idxH'
  where
    (idxG, idxH) = w_black_pair bw
    (idxG', idxH') = w_black_pair bw'

-- return true if the black edges from the witnesses are parallel
checkParallel :: BlockWitness -> BlockWitness -> Bool
checkParallel bw bw' = aligned && intersected
  where
    aligned =
      abs (idxG - idxG') == abs (idxH - idxH')
        && w_inverted bw == w_inverted bw'
    intersected =
      if w_inverted bw
        then
          (idxG' >= idxG && idxH' <= idxH && idxG' <= idxG_r && idxH' >= idxH_r)
            || (idxG' <= idxG && idxH' >= idxH && idxG' >= idxG_l && idxH' <= idxH_l)
        else
          (idxG' >= idxG && idxH' >= idxH && idxG' <= idxG_r && idxH' <= idxH_r)
            || (idxG' <= idxG && idxH' <= idxH && idxG' >= idxG_l && idxH' >= idxH_l)
    (idxG, idxH) = w_black_pair bw
    (idxG_l, idxH_l) = w_left_pair bw
    (idxG_r, idxH_r) = w_right_pair bw
    (idxG', idxH') = w_black_pair bw'

-- check if every vertex has degree at least low and at most up, and it is connect only
-- with vertices of the other genome. In case of black edges, the degree must be
-- exactly 1.
checkDegrees :: (Genome g1, Genome g2) => Int -> Int -> SampleGraph m g1 g2 -> Bool
checkDegrees low up SG {..} =
  all (checkDegree . Vtx True) [1 .. mkIdx (size sg_g)]
    && all (checkDegree . Vtx False) [1 .. mkIdx (size sg_h)]
  where
    checkDegree v =
      let l = edges MMap.! v
          (low', up') = if any eBlack l then (1, 1) else (low, up)
       in low' <= length l
            && length l <= up'
            && not (any ((vInG v ==) . vInG . eOther) l)

edgesFromRareOddPath :: (Matcher m g1 g2, Genome g1, Genome g2) => SampleGraph m g1 g2 -> [Vertex] -> [BlockWitness]
edgesFromRareOddPath sg@(SG {..}) vertices = do
  v <- if length verticesG > length verticesH then verticesG else verticesH
  filter (not . checkBWPairs) . getNewEdges $ v
  where
    checkBWPairs bw = any (\bw' -> checkParallel bw bw' || checkConflict bw bw') sample
    verticesG = filter vInG vertices
    verticesH = filter (not . vInG) vertices
    getNewEdges v =
      map (findLimits sg) $
        if vInG v
          then (vIdx v,) <$> getVertexInQ sg_h (getGene (vIdx v) sg_g)
          else (,vIdx v) <$> getVertexInQ sg_g (getGene (vIdx v) sg_h)
    -- return idx of Q including a vertex with label equal to the one of gene
    getVertexInQ q gene =
      filter
        ( \idx ->
            canonicOri (getGene idx q) == canonicOri gene
        )
        [1 .. mkIdx (size q)]

-- let vsH = map (vIdx . eOther) $ edges MMap.! Vtx True idxG
--  in (idxH `elem` vsH)

createSampleGraph :: (Genome g1, Genome g2) => m g1 g2 -> g1 -> g2 -> Sample -> SampleGraph m g1 g2
createSampleGraph matcher g h = foldr addBlackEdge sg0
  where
    sg0 = SG matcher g h [] MMap.empty

-- Convert pair of indices into a block witness by finding the limits of
-- pairs parallel to it.
findLimits :: (Matcher m g1 g2, Genome g1, Genome g2) => SampleGraph m g1 g2 -> (Idx, Idx) -> BlockWitness
findLimits SG {..} (idxG0, idxH0) =
  BW (idxGL, idxHL) (idxG0, idxH0) (idxGR, idxHR) (not dirH0)
  where
    dirH0 =
      isDirectMatch
        sg_matcher
        (subGenome idxG0 idxG0 sg_g)
        (subGenome idxH0 idxH0 sg_h)
    (idxGL, idxHL) = findLimit False (not dirH0) (idxG0, idxH0)
    (idxGR, idxHR) = findLimit True dirH0 (idxG0, idxH0)

    -- The directions (dir) indicate what is the next gene to look at.
    -- If dir is False we are going from right to left if it is True we are going from
    -- left to right.
    findLimit dirG dirH (idxG, idxH) =
      if a1 < 1
        || a2 > mkIdx (size sg_g)
        || b1 < 1
        || b2 > mkIdx (size sg_h)
        || dirG == dirH
          && not
            (isDirectMatch sg_matcher (subGenome a1 a2 sg_g) (subGenome b1 b2 sg_h))
        || dirG /= dirH
          && not
            (isReverseMatch sg_matcher (subGenome a1 a2 sg_g) (subGenome b1 b2 sg_h))
        then (idxG, idxH)
        else findLimit dirG dirH (idxG', idxH')
      where
        (a1, a2, idxG') =
          if dirG
            then (idxG, incIdx idxG, incIdx idxG)
            else (decIdx idxG, idxG, decIdx idxG)
        (b1, b2, idxH') =
          if dirH
            then (idxH, incIdx idxH, incIdx idxH)
            else (decIdx idxH, idxH, decIdx idxH)

-- Given a correspondence between a gene of G and a gene of H, include in the graph
-- the correspondent black edge and the edges parallel to that correspondence.
-- Also remove any edge that is no longer valid.
addBlackEdge :: (Genome g1, Genome g2) => BlockWitness -> SampleGraph m g1 g2 -> SampleGraph m g1 g2
addBlackEdge bw sg@(SG {..}) = sg {sample = sample', edges = edges'''}
  where
    (idxG0, idxH0) = w_black_pair bw
    sample' = bw : sample
    vG0 = Vtx True idxG0
    vH0 = Vtx False idxH0
    edges_clean = cleanEdges vH0 . cleanEdges vG0 $ edges
    edges' =
      MMap.insert vH0 (Edge vG0 (idxG0, idxH0) bw True) $
        MMap.insert vG0 (Edge vH0 (idxG0, idxH0) bw True) edges_clean
    edges'' = addEdges True (not (w_inverted bw)) (idxG0, idxH0) edges'
    edges''' = addEdges False (w_inverted bw) (idxG0, idxH0) edges''

    -- Remove any green edge that is no longer valid after the inclusion of a new black
    -- edge incident on v.
    cleanEdges v aux_edges = aux_edges'
      where
        aux_edges' = foldr (delGreens v) aux_edges (aux_edges MMap.! v)
    delGreens v edge = delFromVertex (eOther edge) . delFromVertex v
      where
        bw' = eParent edge
        (idxG, idxH) = w_black_pair bw'
        delFromVertex u_aux = go u_aux
          where
            move =
              -- how to move to the next vertex
              if vInG u_aux && vIdx u_aux < idxG || not (vInG u_aux) && vIdx u_aux < idxH
                then \u -> Vtx (vInG u) (decIdx . vIdx $ u)
                else \u -> Vtx (vInG u) (incIdx . vIdx $ u)
            go u aux_edges =
              if vIdx u <= 0
                || (vInG u && vIdx u > mkIdx (size sg_g))
                || (vInG u && vIdx u > mkIdx (size sg_h))
                || length e_list == length e_list'
                then aux_edges
                else go (move u) aux_edges'
              where
                e_list = aux_edges MMap.! u
                e_list' = filter (\e -> eParent e /= bw') e_list
                aux_edges' = MMap.fromMap . Map.insert u e_list' . MMap.toMap $ aux_edges

    -- The directions (dir) indicate what is the next gene to look at.
    -- If dir is False we are going from right to left if it is True we are going from
    -- left to right.
    addEdges dirG dirH (idxG, idxH) aux_edges =
      if not dirG && idxG == fst (w_left_pair bw) || dirG && idxG == fst (w_right_pair bw) || hasBlackEdge vG || hasBlackEdge vH
        then aux_edges
        else addEdges dirG dirH (idxG', idxH') aux_edges'
      where
        hasBlackEdge v = let l = aux_edges MMap.! v in not (null l) && eBlack (head l)
        idxG' = (if dirG then incIdx else decIdx) idxG
        idxH' = (if dirH then incIdx else decIdx) idxH
        vG = Vtx True idxG'
        vH = Vtx False idxH'
        aux_edges' =
          MMap.insert vH (Edge vG (idxG', idxH') bw False) $
            MMap.insert vG (Edge vH (idxG', idxH') bw False) aux_edges

-- Make sure that there is at most one edge incident in each pair of vertices
-- so that we have a correspondence between the rare genes of G and H.
-- At this point the graph should have only black edges and even green paths, so
-- we eliminate green edges to ensure the correspondence.
finalizeCorrespondence :: (Genome g1) => SampleGraph m g1 g2 -> SampleGraph m g1 g2
finalizeCorrespondence sg = sg2
  where
    (seen1, sg1) = foldr (changePaths . Vtx True) (Set.empty, sg) [1 .. (mkIdx . size . sg_g $ sg)]
    sg2 = snd $ foldr (changeCycle . Vtx True) (seen1, sg1) [1 .. (mkIdx . size . sg_g $ sg1)]

    changeCycle v (seen, sg_aux) =
      if v `Set.member` seen
        then (seen, sg_aux)
        else case edges sg_aux MMap.! v of
          [e1, e2] ->
            let e =
                  if (eOther e1 `Set.member` seen)
                    || (fst (w_black_pair . eParent $ e2) < fst (eVpair e2))
                      && (eOther e2 `Set.notMember` seen)
                    then e2
                    else e1
             in changePathsEven (eOther e) (seen, sg_aux)
          _ -> error patternError

    changePaths v (seen, sg_aux) =
      if v `Set.member` seen
        then (seen, sg_aux)
        else case edges sg_aux MMap.! v of
          [_, _] -> (seen, sg_aux) -- not an extreme or is part of a cycle, so we ignore for now
          [e] -> changePathsEven (eOther e) (aux_seen, sg_aux)
          _ -> error patternError
      where
        aux_seen = Set.insert v seen

    changePathsOdd v (seen, sg_aux) =
      if v `Set.member` seen
        then (seen, sg_aux)
        else case edges sg_aux MMap.! v of
          [e] -> changePathsEven (eOther e) (aux_seen, sg_aux)
          _ -> error patternError
      where
        aux_seen = Set.insert v seen

    changePathsEven v (seen, sg_aux) =
      case edges sg_aux MMap.! v of
        [_] -> (aux_seen, sg_aux) -- path has only one vertex remaining
        [e1, e2] ->
          let u =
                if eOther e1 `Set.member` seen
                  then eOther e2
                  else eOther e1
           in changePathsOdd u (aux_seen, delete_edge u v sg_aux)
        _ -> error patternError
      where
        aux_seen = Set.insert v seen
        delete_edge u1 u2 sg' = sg' {edges = updateEdgeInMap u1 u2 . updateEdgeInMap u2 u1 $ edges sg'}
        updateEdgeInMap u1 u2 m =
          MMap.fromMap . Map.insert u1 (delEdge [] (m MMap.! u1)) . MMap.toMap $ m
          where
            delEdge _ [] = error patternError
            delEdge acc (e : es) = if eOther e == u2 then reverse acc ++ es else delEdge (e:acc) es
