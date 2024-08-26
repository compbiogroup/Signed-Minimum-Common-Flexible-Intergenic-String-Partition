{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeFamilies #-}

-- |
-- Module      : DB
-- Description : Code to generate database with artificial genomes represented by
-- strings plus some intergenic region information.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module DB where

import Control.Monad (replicateM, zipWithM_)
import Control.Monad.Random (Rand, StdGen, evalRandIO, getRandom, getRandomR, getRandoms, uniform)
import Data.ByteString (ByteString)
import Data.ByteString.Char8 qualified as BS
import Data.List (transpose, unfoldr)
import Data.Set qualified as Set
import Genomes
import Options.Applicative

opts :: ParserInfo Args
opts =
  info
    (argParser <**> helper)
    ( fullDesc
        <> progDesc "Generate database with genomes. Used for tests of solutions for rearrangement problems."
    )

data Parameters = DB1 RepDB | DB2 RandDB

data RepDB = RepDB
  { db_low :: Int,
    db_high :: Int,
    db_rep :: Int
  }

newtype RandDB = RandDB {db_lim :: Int}

data FlexLevel = Rigid | Flex Int Int

data Args = Args
  { db_par :: Parameters,
    db_flexLevel :: FlexLevel,
    db_num_pairs :: Int,
    db_size :: Int,
    db_nop :: Int,
    db_porc :: Int,
    db_sign :: Sign,
    db_indel :: Int,
    db_mult_chrom :: Bool,
    db_zeros :: Bool,
    db_output :: String
  }

repDBParser :: Parser RepDB
repDBParser =
  RepDB
    <$> option
      auto
      ( long "low"
          <> short 'l'
          <> metavar "L"
          <> help "Minimum number of replicas."
      )
    <*> option
      auto
      ( long "high"
          <> short 'h'
          <> metavar "H"
          <> help "Maximum number of replicas."
      )
    <*> option
      auto
      ( long "replicas"
          <> short 'r'
          <> metavar "REP"
          <> help "Number of genes having more than one occurrence."
      )

randDBParser :: Parser RandDB
randDBParser =
  RandDB
    <$> option
      auto
      ( long "alphabet-size"
          <> short 'l'
          <> metavar "ALPH"
          <> help "Alphabet size."
      )

flexParser :: Parser FlexLevel
flexParser =
  ( \l h ->
      if
          | l == -1 && h == -1 -> Rigid
          | h == -1 -> Flex l l
          | l == -1 -> Flex h h
          | otherwise -> Flex l h
  )
    <$> option
      auto
      ( long "flex-low"
          <> short 'f'
          <> metavar "LF"
          <> value (-1)
          <> help "Whether to produce flexible intergenic regions on the target genome, the value represents the percentage of the intergenic regions size subtracted or added to produce the bound of the interval. If flex-hight is also present generate multiple files with flexibility values flex-low,flex-low + 10,...,flex-hight"
      )
    <*> option
      auto
      ( long "flex-hight"
          <> short 'F'
          <> metavar "HF"
          <> value (-1)
          <> help "Whether to produce flexible intergenic regions on the target genome, the value represents the percentage of the intergenic regions size subtracted or added to produce the bound of the interval. If flex-low is also present generate multiple files with flexibility values flex-low,flex-low + 10,...,flex-hight"
      )

argParser :: Parser Args
argParser =
  Args
    <$> subparser
      ( command
          "RepDB"
          ( info
              (DB1 <$> repDBParser <**> helper)
              (progDesc "DB with fix number of replicas.")
          )
          <> command
            "RandDB"
            ( info
                (DB2 <$> randDBParser <**> helper)
                (progDesc "DB with random generation.")
            )
      )
    <*> flexParser
    <*> option
      auto
      ( long "number_genomes"
          <> short 'k'
          <> metavar "K"
          <> help "Number genome pairs to generate."
      )
    <*> option
      auto
      ( long "size_genome"
          <> short 'n'
          <> metavar "N"
          <> help "Size of the genomes."
      )
    <*> option
      auto
      ( long "number_op"
          <> short 'r'
          <> metavar "R"
          <> help "Number of operations to apply (-1 to use a random list)."
      )
    <*> option
      auto
      ( long "porcentage_rev"
          <> short 'p'
          <> metavar "P"
          <> showDefault
          <> value 100
          <> help "Porcentage of reversions in the operations."
      )
    <*> flag
      Unsigned
      Signed
      ( long "signed"
          <> short 's'
          <> help "Whether the input Strings are signed."
      )
    <*> option
      auto
      ( long "indel"
          <> short 'd'
          <> metavar "D"
          <> help "How much indels to apply (D deletions follow by D insertions)."
      )
    <*> switch
      ( long "mult_chrom"
          <> short 'c'
          <> help "Whether to produce multiple chromosomes after operations. In that case all operations will be DCJs."
      )
    <*> switch
      ( long "zeros"
          <> short 'z'
          <> help "Whether to produce intergenic regions with zeros."
      )
    <*> strOption
      ( long "outfile"
          <> short 'o'
          <> metavar "oFILE"
          <> help "Output file"
      )

main :: IO ()
main = do
  args <- execParser opts
  let (flex_levels, names) = case db_flexLevel args of
        Rigid -> ([Rigid], [db_output args])
        Flex l h ->
          if l == h
            then ([Flex l h], [db_output args])
            else unzip . fmap (\f -> (Flex f f, db_output args ++ "_" ++ show f)) $ [l, l + 10 .. h]
  dbs <- evalRandIO . fmap (map (BS.unlines . fromQuadruples) . transpose) $ replicateM (db_num_pairs args) (genPair args flex_levels)
  zipWithM_ BS.writeFile names dbs
  where
    fromQuadruples ((s1, i1, s2, i2) : ss) = s1 : i1 : s2 : i2 : fromQuadruples ss
    fromQuadruples [] = []

-- TODO: Produce genomes with dcj and indels
genPair :: Args -> [FlexLevel] -> Rand StdGen [(ByteString, ByteString, ByteString, ByteString)]
genPair Args {..} flex_levels = do
  g <- case db_par of
    (DB1 (RepDB rl rh d)) -> randomGenomeWithReplicas db_zeros db_size d rl rh db_sign
    (DB2 (RandDB lim)) -> randomGenome db_zeros db_size lim db_sign
  -- floor of nop for orig and ceil for target
  let nop_orig = db_nop `div` 2
      nop_target = db_nop - nop_orig
  g_final <- if db_mult_chrom && db_nop > 0 then applyDCJs (toMC g) nop_orig else return (toMC g)
  h <-
    if db_mult_chrom && db_nop /= 0
      then (`applyDCJs` nop_target) =<< (toMC <$> applyIndels g)
      else
        toMC
          <$> ( applyIndels
                  =<< if db_nop == -1
                    then shuffleGenome g
                    else applyOperations g
              )
  return $ genBS g_final h
  where
    toMC g' = mkMCRGenome [getGene 1 g', getGene (mkIdx (size g')) g'] [g']
    genBS g h =
      map
        ( \f ->
            let (s2, i2) =
                  ( case f of
                      Rigid -> writeMultiC (writeRGenome True) h
                      Flex l _ -> writeMultiC (writeFGenome True) (flexibilizeMC l h)
                  )
             in (s1, i1, s2, i2)
        )
        flex_levels
      where
        (s1, i1) = writeMultiC (writeRGenome True) g
    
    -- floor of nop for r_r and ceil for r_t
    r_r = (db_nop * db_porc) `div` 100
    r_t = db_nop - r_r

    applyIndels g = do
      let dels = unfoldr dels_for_one db_indel
      g' <- foldr (=<<) (return g) dels
      let ins = unfoldr ins_for_one (db_indel, (+ 1) . geneToInt . (!! 1) . Set.toDescList . alphabet $ g)
      foldr (=<<) (return g') ins
    dels_for_one 0 = Nothing
    dels_for_one d = Just . (,d - 1) $ \g -> do
      i <- getRandomR (2 :: Idx, mkIdx $ size g - 1)
      ir <- getRandomR (0, irToInt (getIR (i - 1) g) + irToInt (getIR i g))
      return $ intergenicDeletion i (i + 1) ir g
    ins_for_one (0, _) = Nothing
    ins_for_one (d, next) = Just . (,(d - 1, next + 1)) $ \g -> do
      i <- getRandomR (1 :: Idx, mkIdx $ size g - 1)
      ir1 <- if db_zeros then return 0 else getRandomR (0, 100)
      ir2 <- if db_zeros then return 0 else getRandomR (max 0 (irToInt (getIR i g) - ir1), 100)
      return $ intergenicInsertion i (mkRGenome Linear False db_sign [next] [ir1, ir2]) g

    applyOperations :: (RigidIntergenicChromosome g) => g -> Rand StdGen g
    applyOperations g = do
      coins <- getRandoms
      let ops = unfoldr operations_for_one (r_t, r_r, coins)
      foldr (=<<) (return g) ops

    applyDCJs :: MultiGIs GenesIRsR -> Int -> Rand StdGen (MultiGIs GenesIRsR)
    applyDCJs g 0 = return g
    applyDCJs g nop = do
      let crhs_idxs = filter (\idx -> let chr = getChromosome idx g in size chr > 1 || isCircular chr && size chr > 0) [1..mkIdx (numChromosomes g)]
      chri_ <- uniform crhs_idxs
      chrj_ <- uniform . filter (\idx -> let chr = getChromosome idx g in idx /= chri_ || size chr > 2 || size chr == 2 && isCircular chr) $ crhs_idxs
      let chri = min chri_ chrj_
          chrj = max chri_ chrj_
          chr_i = getChromosome chri g
          chr_j = getChromosome chrj g
      i <- getRandomR (1, mkIdx $ size chr_i - (if isCircular chr_i then 0 else 1) - (if chri == chrj then 1 else 0))
      j <- getRandomR (if chri == chrj then i+1 else 1, mkIdx $ size chr_j - (if isCircular chr_j then 0 else 1))
      x <- getRandomR (0, irToInt $ getIR i chr_i)
      y <- getRandomR (0, irToInt $ getIR j chr_j)
      connectAC <- getRandom
      applyDCJs (intergenicDCJ (chri, i) (chrj, j) connectAC x y g) (nop - 1)

    operations_for_one :: (RigidIntergenicChromosome g) => (Int, Int, [Bool]) -> Maybe (g -> Rand StdGen g, (Int, Int, [Bool]))
    operations_for_one (_, _, []) = Nothing
    operations_for_one (r_t', r_r', coin : coins)
      | r_t' == 0 && r_r' == 0 = Nothing
      | r_t' == 0 || r_r' /= 0 && coin = Just . (,(r_t', r_r' - 1, coins)) $ \g -> do
          i_ <- getRandomR (2, db_size - 2)
          j_ <- getRandomR (i_ + 1, db_size - 1)
          let (i, j) = (mkIdx i_, mkIdx j_)
          x <- getRandomR (0, irToInt $ getIR (i - 1) g)
          y <- getRandomR (0, irToInt $ getIR j g)
          return $ intergenicReversal i j x y g
      | otherwise = Just . (,(r_t' - 1, r_r', coins)) $ \g -> do
          i_ <- getRandomR (2, db_size - 2)
          j_ <- getRandomR (i_ + 1, db_size - 1)
          k_ <- getRandomR (j_ + 1, db_size)
          let (i, j, k) = (mkIdx i_, mkIdx j_, mkIdx k_)
          x <- getRandomR (0, irToInt $ getIR (i - 1) g)
          y <- getRandomR (0, irToInt $ getIR (j - 1) g)
          z <- getRandomR (0, irToInt $ getIR (k - 1) g)
          return $ intergenicTransposition i j k x y z g
