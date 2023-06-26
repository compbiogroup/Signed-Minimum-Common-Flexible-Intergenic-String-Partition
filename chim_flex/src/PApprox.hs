{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE TupleSections #-}

module PApprox
  ( approxPartition,
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
import Data.IntSet qualified as IntSet
import Data.Map (Map)
import Data.Map qualified as Map
import Data.Maybe (isJust, listToMaybe, mapMaybe)
import Data.Set (Set)
import Data.Set qualified as Set
import Genomes (Gene, Genome (getGene, size, subGenome), Idx, Matcher (isMatch), mkIdx)
import LocalBase
import Partition (CommonPartition, Partition, blocks, breakpoints, mkCommonPartition, mkPartition, trivialPartition, underlineGenome)

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

data BlockMatchGraph g1 g2 = BG (Partition g1) (Partition g2) (Map Vertex (Maybe Int, Correspondents)) deriving (Show)

data Component = Component
  { cVtx :: Vertex,
    cContG :: Int,
    cContH :: Int,
    cPerfect :: Bool
  }
  deriving (Show)

countPerfect :: IntMap Component -> Int
countPerfect = sum . map (fromEnum . cPerfect) . IntMap.elems

approxPartition :: (Orientable g1, Orientable g2, Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> g1 -> g2 -> CommonPartition g1 g2
approxPartition matcher g h = mkCommonPartition pg ph
  where
    pg_0 = trivialPartition g
    ph_0 = trivialPartition h
    (comps, bmg) = getConnectedComponents (mkBlockMatchGraph matcher pg_0 ph_0)
    breakSet_0 = checkIntersections comps bmg

    (_, BG pg ph _, _, _) = until done (addBreaks matcher) (comps, bmg, breakSet_0, False)
    done (_, _, _, b) = b

addBreaks :: (Genome g1, Genome g2, Matcher m g1 g2) => m g1 g2 -> (IntMap Component, BlockMatchGraph g1 g2, Set (Gene, Gene), Bool) -> (IntMap Component, BlockMatchGraph g1 g2, Set (Gene, Gene), Bool)
addBreaks matcher (comps, bmg@(BG pg ph corres), breakSet, _) =
  if null notPerfs then (comps, bmg, breakSet, True) else (comps', bmg', breakSet, False)
  where
    -- The substrings that may recieve a breakpoint are smallest in components that do not admit a perfect match (this ensures that they are in Tmin)
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
    pg' = mkPartition g bg'
    ph' = mkPartition h bh'
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
      let (a, b) = (canonicOri *** canonicOri) $
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
              perfect = isJust (testPerfectMatch corres''' vs)
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
    -- checks if c admits a perfect match
    testPerfectMatch aux_corres = foldrM (selectPBlock Set.empty) Map.empty
      where
        -- blocksFromG saves blocks from G correspondent to each position of blocks in H
        selectPBlock seen v blocksFromG =
          if done then Just blocksFromG' else Nothing
          where
            (_, blocksFromG', done) = selectPBlock' seen v blocksFromG

        selectPBlock' seen v = testCandidates v neigs seen
          where
            (_, neigs) = aux_corres Map.! v

        testCandidates _ [] seen blocksFromG = (seen, blocksFromG, False)
        testCandidates v (u : us) seen blocksFromG =
          if
              | u `Set.member` seen -> testCandidates v us seen blocksFromG
              | u `Map.notMember` blocksFromG -> (seen', blocksFromG', True)
              | done -> (seen_rec, blocksFromG_rec', True)
              | otherwise -> testCandidates v us seen_rec blocksFromG_rec
          where
            seen' = Set.insert u seen
            blocksFromG' = Map.insert u v blocksFromG
            (seen_rec, blocksFromG_rec, done) = selectPBlock' seen' (blocksFromG Map.! u) blocksFromG
            blocksFromG_rec' = Map.insert u v blocksFromG_rec

-- Construct the block match graph, which is represented by an adjacency list decorated -- with information about the component (to be filled by getConnectedComponents later).
-- As we are considering balanced genomes, we ignore subgenomes of size 1.
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

checkIntersections :: (Orientable g1, Orientable g2, Genome g1, Genome g2) => IntMap Component -> BlockMatchGraph g1 g2 -> Set (Gene, Gene)
checkIntersections comps (BG pg ph corres) = do
  interMapToBreakSet . foldr updateInterMap interMap_0 $ vertexPairsTmin
  where
    g = underlineGenome pg
    h = underlineGenome ph

    -- Vertices from components that do not admit a perfect match
    notPerfs = filter (\(_, c) -> not . cPerfect $ comps IntMap.! c) . map getComp . Map.toAscList $ corres
    vertexPairsNotPerfs = [(v_cv, u_cu) | v_cv <- notPerfs, u_cu <- notPerfs]
    getComp (_, (Nothing, _)) = error patternError
    getComp (v, (Just c, _)) = (v, c)
    -- Vertices from components of Tmin
    tmin = filter (\(_, c) -> c `IntSet.member` tminSet) notPerfs
    -- We ensure in these pairs that v starts before u
    vertexPairsTmin = [((v, cv), (u, cu)) | (v, cv) <- tmin, (u, cu) <- tmin, vIdx v < vIdx u]
    tminSet = foldr selectTmin (IntSet.fromList (map snd notPerfs)) vertexPairsNotPerfs

    interMap_0 :: IntMap (Maybe Idx, Maybe Idx)
    interMap_0 = IntMap.fromAscList . map (,(Nothing, Nothing)) . IntSet.toAscList $ tminSet

    selectTmin ((v, cv), (u, cu)) tminSet_aux =
      if
          | cv == cu || vInG v /= vInG u -> tminSet_aux
          | (vBeg v <= vBeg u && vEnd u <= vEnd v) -> IntSet.delete cv tminSet_aux
          | (vBeg u <= vBeg v && vEnd v <= vEnd u) -> IntSet.delete cu tminSet_aux
          | otherwise -> tminSet_aux

    updateInterMap :: ((Vertex, Int), (Vertex, Int)) -> IntMap (Maybe Idx, Maybe Idx) -> IntMap (Maybe Idx, Maybe Idx)
    updateInterMap ((v, cv), (u, cu)) interMap =
      if cv == cu || vInG v /= vInG u || vEnd v <= vBeg u
        then interMap -- Ignore cases where both vertices are in the same component or if there is no intersection
        else
          let interMap_up1 = updateEntry cu (vOri u) lu ru lu (Just $ vEnd u - vEnd v + 1) interMap
           in updateEntry cv (vOri v) lv rv (Just $ vBeg u - vBeg v + 1) rv interMap_up1
      where
        -- For each component, the maps stores the smallest size of the intersection
        -- between components on the left (lv) and on the right (rv)
        (lv, rv) = interMap IntMap.! cv
        (lu, ru) = interMap IntMap.! cu
        vOri v' =
          if vInG v'
            then getOri (subGenome (vBeg v') (vEnd v') g)
            else getOri (subGenome (vBeg v') (vEnd v') h)

        -- update the intersection map taking into account that blocks of a component -- may be inverted
        updateEntry c' ori l_old r_old l_new r_new interMap' =
          if ori == vOri v'
            then IntMap.insert c' (select l_old l_new, select r_old r_new) interMap'
            else IntMap.insert c' (select l_old r_new, select r_old l_new) interMap'
          where
            v' = cVtx $ comps IntMap.! c'
            select old new = min <$> old <*> new <|> old <|> new

    interMapToBreakSet :: IntMap (Maybe Idx, Maybe Idx) -> Set (Gene, Gene)
    interMapToBreakSet = foldr chooseBreak Set.empty . IntMap.toList

    chooseBreak :: (Int, (Maybe Idx, Maybe Idx)) -> Set (Gene, Gene) -> Set (Gene, Gene)
    chooseBreak (c, lr) breakSet =
      case lr of
        (Nothing, Nothing) -> Set.insert (getGenePair (vBeg v)) breakSet
        (Just l, Nothing) -> Set.insert (getGenePair (vBeg v + l - 2)) breakSet
        (Nothing, Just r) -> Set.insert (getGenePair (vEnd v - r + 1)) breakSet
        (Just l, Just r) ->
          let l_inter = vBeg v + l - 1
              r_inter = vEnd v - r + 1
           in if l_inter > r_inter
                then Set.insert (getGenePair r_inter) breakSet
                else
                  Set.insert
                    (getGenePair $ vBeg v)
                    (Set.insert (getGenePair (vEnd v - 1)) breakSet)
      where
        v = cVtx $ comps IntMap.! c
        getGenePair i =
          canonicOri *** canonicOri $
            if vInG v
              then (getGene i g, getGene (i + 1) g)
              else (getGene i h, getGene (i + 1) h)