{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TupleSections #-}

-- |
-- Module      : DB
-- Description : Code to generate database with artificial genomes represented by
-- strings plus some intergenic region information.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module DB where

import Control.Monad (replicateM)
import Control.Monad.Random (Rand, StdGen, evalRandIO, getRandomR, getRandoms)
import Data.ByteString (ByteString)
import Data.ByteString.Char8 qualified as BS
import Data.List (unfoldr)
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

data FlexLevel = Rigid | Flex Int

data Args = Args
  { db_par :: Parameters,
    db_flexLevel :: FlexLevel,
    db_num_pairs :: Int,
    db_size :: Int,
    db_nop :: Int,
    db_porc :: Int,
    db_sign :: Sign,
    db_indel :: Int,
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
  (\i -> if i == -1 then Rigid else Flex i)
    <$> option
      auto
      ( long "flex"
          <> short 'f'
          <> metavar "f"
          <> value (-1)
          <> help "Whether to produce flexible intergenic regions on the target genome, the value represents the percentage of the intergenic regions size subtracted or added to produce the bound of the interval."
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
  db <- evalRandIO . fmap (BS.unlines . fromQuadruples) $ replicateM (db_num_pairs args) (genPair args)
  BS.writeFile (db_output args) db
  where
    fromQuadruples ((s1, i1, s2, i2) : ss) = s1 : i1 : s2 : i2 : fromQuadruples ss
    fromQuadruples [] = []

genPair :: Args -> Rand StdGen (ByteString, ByteString, ByteString, ByteString)
genPair Args {..} = do
  g <- case db_par of
    (DB1 (RepDB l h d)) -> randomGenomeWithReplicas db_zeros db_size d l h db_sign
    (DB2 (RandDB lim)) -> randomGenome db_zeros db_size lim db_sign
  h <-
    applyIndels
      =<< if db_nop == -1
        then shuffleGenome g
        else applyOperations g
  let (s1, i1) = writeRGenome True g
  let (s2, i2) = case db_flexLevel of
        Rigid -> writeRGenome True h
        Flex l -> writeFGenome True (flexibilize l h)
  return (s1, i1, s2, i2)
  where
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
      ir2 <- if db_zeros then return 0 else getRandomR (max 0 ((irToInt $ getIR i g) - ir1), 100)
      return $ intergenicInsertion i (mkRGenome False db_sign [next] [ir1, ir2]) g

    applyOperations :: (RigidIntergenicGenome g, IntergenicGenome g) => g -> Rand StdGen g
    applyOperations g = do
      coins <- getRandoms
      let ops = unfoldr operations_for_one (r_t, r_r, coins)
      foldr (=<<) (return g) ops

    operations_for_one :: (RigidIntergenicGenome g, IntergenicGenome g) => (Int, Int, [Bool]) -> Maybe (g -> Rand StdGen g, (Int, Int, [Bool]))
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
