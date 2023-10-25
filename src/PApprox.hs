{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE TupleSections #-}

-- \|
-- Module      : PApprox
-- Description : Approximation algorithm for partition problems. The approximation
-- factor depends on the maximum occurrence of a gene in the genomes.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module PApprox
  ( approxPartition,
    approxLowerBound,
    countPerfect,
    getConnectedComponents,
    mkBlockMatchGraph,
  )
where

import Control.Applicative ((<|>))
import Control.Arrow ((***))
import Data.EnumSet qualified as EnumSet
import Data.Foldable (foldrM)
import Data.IntMap (IntMap)
import Data.IntMap qualified as IntMap
import Data.IntSet (IntSet)
import Data.IntSet qualified as IntSet
import Data.Map (Map)
import Data.Map qualified as Map
import Data.Maybe (isJust, listToMaybe, mapMaybe)
import Data.Set (Set)
import Data.Set qualified as Set
import Genomes (FlipMatcher (FlipMatcher), Gene, Genome (getGene, size, subGenome), Idx, Matcher (isDirectMatch, isMatch), isCompatibleWithSubGenome, mkIdx)
import LocalBase
import Partition (CommonPartition, Partition, blocks, breakpoints, findUncommon, mkCommonPartition, mkPartitionFromBlocks, mkPartitionFromBreakpoints, trivialPartition, underlineGenome)

type Correspondents = [Vertex]

data Vertex = Vtx
  { vSize :: Int,
    vInG :: Bool,
    vIdx :: Idx
  }
  deriving (Eq, Ord, Show)

vBeg :: Vertex -> Idx
vBeg = vIdx

vEnd :: Vertex -> Idx
vEnd v = vIdx v + mkIdx (vSize v) - 1

-- Graph where each vertex is a block and an edge connects two vertices if they are correspondents, also has a space to record the number of the connected components containing each graph
data BlockMatchGraph g1 g2 = BG (Partition g1) (Partition g2) (Map Vertex (Maybe Int, Correspondents)) deriving (Show)

data Component = Component
  { cVtx :: Vertex,
    cContG :: Int,
    cContH :: Int,
    cPerfect :: Bool
  }
  deriving (Show)

-- count the number of perfect components
countPerfect :: IntMap Component -> Int
countPerfect = sum . map (fromEnum . cPerfect) . IntMap.elems

approxPartition :: (Orientable g1, Orientable g2, Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> g1 -> g2 -> CommonPartition g1 g2
approxPartition matcher g h = mkCommonPartition matcher pg'' ph''
  where
    pg_0 = trivialPartition g
    ph_0 = trivialPartition h

    (pg'', ph'', _, _) = until done2 tryToBuildPartition (pg_0, ph_0, Nothing, False)
    tryToBuildPartition (pg, ph, maybe_breakSet, _) = compatibilyCompletion matcher pg' ph' breakSet
      where
        (_, BG pg' ph' _) = until done1 (addBreaks matcher breakSet) (comps, bmg)
        done1 (comps', _) = countPerfect comps' == IntMap.size comps'
        (comps, bmg) = getConnectedComponents (mkBlockMatchGraph matcher pg ph)
        breakSet = case maybe_breakSet of
          Just bs -> bs
          Nothing -> fst $ checkIntersections comps bmg
    done2 (_, _, _, isCommon) = isCommon

approxLowerBound :: (Orientable g1, Orientable g2, Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> g1 -> g2 -> Int
approxLowerBound matcher g h = lb
  where
    pg_0 = trivialPartition g
    ph_0 = trivialPartition h
    (comps, bmg) = getConnectedComponents (mkBlockMatchGraph matcher pg_0 ph_0)
    (_, lb) = checkIntersections comps bmg

-- Add breakpoints in some subgenome from Tmin
addBreaks :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> Set (Gene, Gene) -> (IntMap Component, BlockMatchGraph g1 g2) -> (IntMap Component, BlockMatchGraph g1 g2)
addBreaks matcher breakSet (comps, BG pg ph corres) = (comps', bmg')
  where
    -- The substrings that may receive a breakpoint are smallest in components that do not admit a perfect match (this ensures that they are in Tmin)
    cand = takeWhile (\v' -> vSize (fst . head $ notPerfs) == vSize (fst v')) notPerfs
    notPerfs = filter isNotPerfect . Map.toAscList $ corres
    isNotPerfect (_, (Nothing, _)) = error patternError
    isNotPerfect (_, (Just c, _)) = not . cPerfect $ comps IntMap.! c
    v = fst . head . sortWith sortKey $ cand

    -- We give priority to cut a substring in the genome with more vertices of the
    -- component and break ties with the degree of the vertex.
    sortKey (v_, (Just c_id, neigs)) = (posWeight, dg)
      where
        c = comps IntMap.! c_id
        posWeight =
          if
              | cContG c == cContH c -> (1 :: Int)
              | (cContG c > cContH c) == vInG v_ -> 0
              | otherwise -> 2
        dg = length neigs
    sortKey (_, (_, _)) = error patternError

    g = underlineGenome pg
    h = underlineGenome ph
    (comps', bmg') = getConnectedComponents (mkBlockMatchGraph matcher pg' ph')
    bg = breakpoints pg
    bh = breakpoints ph
    pg' = mkPartitionFromBreakpoints g bg'
    ph' = mkPartitionFromBreakpoints h bh'
    (bg', bh') =
      if vInG v
        then (EnumSet.insert break_idx bg, bh)
        else (bg, EnumSet.insert break_idx bh)
    break_idx =
      case testIRs of
        Just b -> b
        Nothing -> error patternError
    testIRs = listToMaybe . mapMaybe testIR $ [vBeg v .. vEnd v - 1]
    testIR i =
      let (a, b) =
            (canonicOri *** canonicOri) $
              if vInG v
                then (getGene i g, getGene (i + 1) g)
                else (getGene i h, getGene (i + 1) h)
       in if (a, b) `Set.member` breakSet || (b, a) `Set.member` breakSet
            then Just i
            else Nothing

-- | Get a set with all connected componets of the Block Match Graph and
-- record the component indices in the Graph
getConnectedComponents :: BlockMatchGraph g1 g2 -> (IntMap Component, BlockMatchGraph g1 g2)
getConnectedComponents (BG pg ph corres) = (comps, BG pg ph corres')
  where
    vertices = Map.keys corres
    (comps, corres', _) = foldr checkComponent (IntMap.empty, corres, 0) vertices

    checkComponent v (comps', corres'', comp_id) =
      case fst (corres'' Map.! v) of
        Just _ -> (comps', corres'', comp_id)
        Nothing ->
          let (corres''', vs, cont_G, cont_H) = visitComp comp_id v (corres'', [], 0, 0)
              perfect = testPerfectMatch corres''' vs cont_G cont_H
              comp = Component v cont_G cont_H perfect
              comps'' = IntMap.insert comp_id comp comps'
           in (comps'', corres''', comp_id + 1)

    visitComp comp_id v (aux_corres, vs, cont_G, cont_H) =
      let (maybe_comp, neigs) = aux_corres Map.! v
       in case maybe_comp of
            Just _ -> (aux_corres, vs, cont_G, cont_H)
            Nothing ->
              let aux_corres' = Map.insert v (Just comp_id, neigs) aux_corres
                  (cont_G', cont_H') =
                    if vInG v
                      then (cont_G + 1, cont_H)
                      else (cont_G, cont_H + 1)
               in foldr (visitComp comp_id) (aux_corres', v : vs, cont_G', cont_H') neigs

    -- Given a correspondence and a list of vertices from a component c,
    -- checks if c admits a perfect match.
    -- Note that we check both directions if the blocks in G can be assigned to blocks
    -- of H and if blocks of H can be assigned to blocks of G.
    testPerfectMatch :: Map Vertex (Maybe Int, Correspondents) -> [Vertex] -> Int -> Int -> Bool
    testPerfectMatch aux_corres allVertices cont_G cont_H =
      cont_G == cont_H && not (null hVertices) && isJust (foldrM (selectBlock Set.empty) Map.empty hVertices)
      where
        -- select only vertices from H, to avoid double checking
        hVertices = filter (not . vInG) allVertices
        -- Select block from G that corresponds to block represented by v.
        -- fixCorres saves a map between blocks of G and their fix correspondent
        -- in H.
        selectBlock seen v fixCorres =
          if success then Just fixCorres' else Nothing
          where
            (_, fixCorres', success) = selectBlock' seen v fixCorres

        selectBlock' seen v = testCandidates v neigs seen
          where
            (_, neigs) = aux_corres Map.! v

        testCandidates _ [] seen fixCorres = (seen, fixCorres, False)
        testCandidates v (u : us) seen fixCorres =
          if
              | u `Set.member` seen -> testCandidates v us seen fixCorres
              | u `Map.notMember` fixCorres -> (seen', fixCorres', True)
              | success -> (seen_rec, fixCorres_rec', True)
              | otherwise -> testCandidates v us seen_rec fixCorres_rec
          where
            seen' = Set.insert u seen
            fixCorres' = Map.insert u v fixCorres
            (seen_rec, fixCorres_rec, success) = selectBlock' seen' (fixCorres Map.! u) fixCorres
            fixCorres_rec' = Map.insert u v fixCorres_rec

-- | Construct the block match graph, which is represented by an adjacency list decorated -- with information about the component (to be filled by getConnectedComponents later).
-- As we are considering balanced genomes, we ignore subgenomes of size 1.
mkBlockMatchGraph :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> Partition g1 -> Partition g2 -> BlockMatchGraph g1 g2
mkBlockMatchGraph matcher pg ph = BG pg ph corres
  where
    bsG = blocks pg
    bsH = blocks ph
    corres = Map.fromList $ corresList bsG getCorresInH True ++ corresList bsH getCorresInG False
    getCorresInH = getCorresIn bsH (isMatch matcher) False
    getCorresInG = getCorresIn bsG (flip (isMatch matcher)) True

    corresList bsK getCorresInK inG = do
      (b, prev_sizes) <- acc_block_sizes bsK
      beg <- [1 .. mkIdx (size b) - 1]
      end <- [beg + 1 .. mkIdx (size b)]
      let sub = subGenome beg end b
      return (Vtx (size sub) inG (mkIdx prev_sizes + beg), (Nothing, getCorresInK sub))

    getCorresIn bsK test inG sub = do
      (b, prev_sizes) <- acc_block_sizes bsK
      beg <- [1 .. mkIdx (size b) - 1]
      let end = beg + mkIdx (size sub) - 1
      [ Vtx (size sub) inG (mkIdx prev_sizes + beg)
        | end <= mkIdx (size b) && test sub (subGenome beg end b)
        ]

    acc_block_sizes :: (Genome g) => [g] -> [(g, Int)]
    acc_block_sizes bsK = zip bsK . scanl (+) 0 . map size $ bsK

-- | Given a set of components and a block match graph, returns a set of breakpoints to be included on the genomes to produce the partitions and
-- a lower bound on the number of breakpoints that are really necessary in a minimum common partition between the genomes.
checkIntersections :: (Orientable g1, Orientable g2, Genome g1, Genome g2) => IntMap Component -> BlockMatchGraph g1 g2 -> (Set (Gene, Gene), Int)
checkIntersections comps (BG pg ph corres) = (Set.unions [breakSet_i, breakSet_ii, breakSet_iii], Set.size breakSet_i + Set.size breakSet_ii + Set.size breakSet_iii `div` 2)
  where
    g = underlineGenome pg
    h = underlineGenome ph

    -- case i: components of tmin that do not have intersections with the others
    -- case ii: components of tmin that with at least of ir that intersect all intersection components
    -- case iii: other components of tmin, must intersect all intersection components in one of two breakpoints (extremes of the subgenomes)
    (cases_i, cases_ii, cases_iii) = foldr separetInterCases ([], [], []) $ IntSet.toList tminSet
    breakSet_i = foldr addBreaks_i Set.empty cases_i
    (breakSet_ii, remaining_cases_iii) = snd $ until (IntSet.null . fst) addBreaks_ii (IntSet.fromList cases_ii, (Set.empty, IntSet.fromList cases_iii))
    breakSet_iii = snd $ until (IntSet.null . fst) addBreaks_iii (remaining_cases_iii, Set.empty)

    -- Vertices from components that do not admit a perfect match
    notPerfs = filter (\(_, c) -> not . cPerfect $ comps IntMap.! c) . map getComp . Map.toAscList $ corres
    vertexPairsNotPerfs = [(v_cv, u_cu) | v_cv <- notPerfs, u_cu <- notPerfs]
    getComp (_, (Nothing, _)) = error patternError
    getComp (v, (Just c, _)) = (v, c)

    -- Vertices from components of Tmin (minimal components without perfect matching)
    tminVtxs = filter (\(_, c) -> c `IntSet.member` tminSet) notPerfs
    -- We ensure in these pairs that v starts before u
    vertexPairsTmin = [((v, cv), (u, cu)) | (v, cv) <- tminVtxs, (u, cu) <- tminVtxs, vIdx v < vIdx u]
    tminSet = foldr selectTmin (IntSet.fromList (map snd notPerfs)) vertexPairsNotPerfs

    -- Starting from a set with all vertices from components without perfect matching,
    -- delete the vertices with components containing subcomponents without perfect matching
    selectTmin ((v, cv), (u, cu)) tminSet_aux =
      if
          | cv == cu || vInG v /= vInG u -> tminSet_aux
          | (vBeg v <= vBeg u && vEnd u <= vEnd v) -> IntSet.delete cv tminSet_aux
          | (vBeg u <= vBeg v && vEnd v <= vEnd u) -> IntSet.delete cu tminSet_aux
          | otherwise -> tminSet_aux

    -- For each component, the intersection Map stores the smallest size of the intersection
    -- between components on the left (lv) and on the right (rv), and the set of components
    -- that intersect it.
    interMap :: IntMap ((Maybe Idx, Maybe Idx), IntSet)
    interMap = foldr updateInterMap interMap_0 vertexPairsTmin
    interMap_0 = IntMap.fromAscList . map (,((Nothing, Nothing), IntSet.empty)) . IntSet.toAscList $ tminSet

    updateInterMap :: ((Vertex, Int), (Vertex, Int)) -> IntMap ((Maybe Idx, Maybe Idx), IntSet) -> IntMap ((Maybe Idx, Maybe Idx), IntSet)
    updateInterMap ((v, cv), (u, cu)) interMap_aux =
      if cv == cu || vInG v /= vInG u || vEnd v <= vBeg u
        then -- Ignore cases where both vertices are in the same component or if there is no intersection
        -- (Note that we already have vIdx v < vIdx u from the creation of vertexPairsTmin)
          interMap_aux
        else
          let interMap_up1 = updateEntry cu cv (vOri u) lu ru interSize ru interMap_aux
           in updateEntry cv cu (vOri v) lv rv lv interSize interMap_up1
      where
        interSize = Just $ vEnd v - vBeg u + 1
        (lv, rv) = fst $ interMap_aux IntMap.! cv
        (lu, ru) = fst $ interMap_aux IntMap.! cu
        vOri v' =
          if vInG v'
            then getOri (subGenome (vBeg v') (vEnd v') g)
            else getOri (subGenome (vBeg v') (vEnd v') h)

        -- update the intersection map taking into account that blocks of a component -- may be inverted
        updateEntry c' c'' ori l_old r_old l_new r_new interMap_aux' =
          if ori == vOri v'
            then IntMap.insert c' ((select l_old l_new, select r_old r_new), IntSet.insert c'' inter) interMap_aux'
            else IntMap.insert c' ((select l_old r_new, select r_old l_new), IntSet.insert c'' inter) interMap_aux'
          where
            v' = cVtx $ comps IntMap.! c'
            select old new = min <$> old <*> new <|> old <|> new
            inter = snd $ interMap_aux' IntMap.! c'

    -- separate connected components of tmin into the three intersection cases
    separetInterCases :: Int -> ([Int], [Int], [Int]) -> ([Int], [Int], [Int])
    separetInterCases c (cs1, cs2, cs3) =
      case fst (interMap IntMap.! c) of
        (Nothing, Nothing) -> (c : cs1, cs2, cs3)
        (Just _, Nothing) -> (cs1, c : cs2, cs3)
        (Nothing, Just _) -> (cs1, c : cs2, cs3)
        (Just l, Just r) ->
          let v = cVtx $ comps IntMap.! c
              l_inter = vBeg v + l - 1
              r_inter = vEnd v - r + 1
           in if l_inter > r_inter
                then (cs1, c : cs2, cs3)
                else (cs1, cs2, c : cs3)

    -- In case i, we can add the breakpoint in any position of the subgenomes
    addBreaks_i :: Int -> Set (Gene, Gene) -> Set (Gene, Gene)
    addBreaks_i c = Set.insert (getGenePair (vBeg v) (vInG v))
      where
        v = cVtx $ comps IntMap.! c

    -- In case ii, for a component c, we add a breakpoint in the intersection off all component intersecting c and update the set of remaining components of tmin
    addBreaks_ii :: (IntSet, (Set (Gene, Gene), IntSet)) -> (IntSet, (Set (Gene, Gene), IntSet))
    addBreaks_ii (cs2, (breakSet, cs3)) = (cs2', (breakSet', cs3'))
      where
        breakSet' = case fst (interMap IntMap.! c) of
          (Nothing, Nothing) -> error patternError
          (Just l, Nothing) -> Set.insert (getGenePair (vBeg v + l - 2) (vInG v)) breakSet
          (_, Just r) -> Set.insert (getGenePair (vEnd v - r + 1) (vInG v)) breakSet
        v = cVtx $ comps IntMap.! c
        c = head $ IntSet.elems cs2
        cs2' = IntSet.delete c cs2 IntSet.\\ snd (interMap IntMap.! c)
        cs3' = IntSet.delete c cs3 IntSet.\\ snd (interMap IntMap.! c)

    -- In case iii, for a component c, we add two breakpoints, such that each component intersecting c will have one of them.
    addBreaks_iii :: (IntSet, Set (Gene, Gene)) -> (IntSet, Set (Gene, Gene))
    addBreaks_iii (cs, breakSet) = (cs', breakSet')
      where
        breakSet' = Set.insert (getGenePair (vBeg v) (vInG v)) (Set.insert (getGenePair (vEnd v - 1) (vInG v)) breakSet)
        v = cVtx $ comps IntMap.! c
        c = head $ IntSet.elems cs
        cs' = IntSet.delete c cs IntSet.\\ snd (interMap IntMap.! c)

    getGenePair i inG =
      canonicOri *** canonicOri $
        if inG
          then (getGene i g, getGene (i + 1) g)
          else (getGene i h, getGene (i + 1) h)

-- | Given two partitions @pg@ and @ph@, with all components admitting a perfect matching,
-- find unmatched blocks and add breakpoints to produce new matches. This function may return
-- partitions still incomplete if some of the matchings are no longer perfect (the boolean in
-- the result is False in this case).
compatibilyCompletion :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> Partition g1 -> Partition g2 -> Set (Gene, Gene) -> (Partition g1, Partition g2, Maybe (Set (Gene, Gene)), Bool)
compatibilyCompletion matcher pg_0 ph_0 breakSet_0 = breakBlock (pg_0, ph_0, breakSet_0)
  where
    -- Find blocks without matches and black bigger blocks to free their matches
    breakBlock (pg, ph, breakSet) =
      case (findUncommon matcher pg ph, findUncommon (FlipMatcher matcher) ph pg) of
        (Just uncommonG, _) ->
          maybe
            (pg, ph, Just breakSet, False)
            breakBlock
            (findBreakpoint matcher breakSet uncommonG pg ph)
        (Nothing, Just uncommonH) ->
          case findBreakpoint (FlipMatcher matcher) breakSet uncommonH ph pg of
            Nothing -> (pg, ph, Just breakSet, False)
            Just (a, b, c) -> breakBlock (b, a, c)
        (Nothing, Nothing) -> (pg, ph, Just breakSet_0, True)

    findBreakpoint matcher' breakSet uncommon p1 p2 =
      case findBreakpoint' uncommon b2s_0 [] of
        Nothing -> Nothing
        Just (breakInfo, new_b2s, new_breaks) ->
          let new_b1s = findBreakpoint'' breakInfo b1s_0 []
              p1' = mkPartitionFromBlocks (underlineGenome p1) new_b1s
              p2' = mkPartitionFromBlocks (underlineGenome p2) new_b2s
              breakSet' = Set.union breakSet (Set.fromList new_breaks)
           in Just (p1', p2', breakSet')
      where
        b1s_0 = blocks p1
        b2s_0 = blocks p2

        -- find block of p2 that contains a substring compatible with a block of uncommon
        findBreakpoint' _ [] _ = Nothing
        findBreakpoint' [] (b2 : b2s) b2s' = findBreakpoint' uncommon b2s (b2 : b2s')
        findBreakpoint' (u : us) (b2 : b2s) b2s' =
          if size u >= size b2
            then findBreakpoint' uncommon b2s (b2 : b2s')
            else case isCompatibleWithSubGenome matcher' u b2 of
              Nothing -> findBreakpoint' us (b2 : b2s) b2s'
              Just i ->
                let b' = subGenome i (i + mkIdx (size u) - 1) b2
                    getGenePair i_ = canonicOri *** canonicOri $ (getGene i_ b2, getGene (i_ + 1) b2)
                    (new_bs, new_breaks) =
                      if
                          | i == 1 -> ([b', subGenome (mkIdx (size u) + 1) (mkIdx (size b2)) b2], [getGenePair (mkIdx (size u))])
                          | i + mkIdx (size u) - 1 == mkIdx (size b2) -> ([subGenome 1 (i - 1) b2, b'], [getGenePair (i - 1)])
                          | otherwise -> ([subGenome 1 (i - 1) b2, b', subGenome (i + mkIdx (size u)) (mkIdx (size b2)) b2], [getGenePair (i - 1), getGenePair (i + mkIdx (size u) - 1)])
                 in Just ((b2, b', i), reverse b2s' ++ new_bs ++ b2s, new_breaks)

        -- find block of p1 that is compatible with the broken block b2 from p2
        findBreakpoint'' (old_b2, new_b, bp_idx) = findBreakpoint''_
          where
            findBreakpoint''_ [] _ = error logicError
            findBreakpoint''_ (b1 : b1s) b1s' =
              if size b1 /= size old_b2 || not (isMatch matcher' b1 old_b2)
                then findBreakpoint''_ b1s (b1 : b1s')
                else
                  let i = if isDirectMatch matcher' b1 old_b2 then bp_idx else mkIdx (size b1) - bp_idx + 1 - mkIdx (size new_b) + 1
                   in let b' = subGenome i (i + mkIdx (size new_b) - 1) b1
                          new_bs =
                            if
                                | i == 1 -> [b', subGenome (mkIdx (size new_b) + 1) (mkIdx (size b1)) b1]
                                | i + mkIdx (size new_b) - 1 == mkIdx (size b1) -> [subGenome 1 (i - 1) b1, b']
                                | otherwise -> [subGenome 1 (i - 1) b1, b', subGenome (i + mkIdx (size new_b)) (mkIdx (size b1)) b1]
                       in reverse b1s' ++ new_bs ++ b1s
