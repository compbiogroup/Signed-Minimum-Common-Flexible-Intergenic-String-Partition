{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module LocalBase where

import qualified Data.HashSet as HashSet
import qualified Data.List as List
import Data.Hashable (Hashable)
import Debug.Trace (trace)

newtype Size = Size Int deriving newtype (Eq, Show, Read, Num, Integral, Real, Ord, Enum)

newtype Dist = Dist Int deriving newtype (Eq, Show, Read)

data Ori = LR | RL deriving (Eq, Show)

class Orientable o where
  -- ^ Get orientation

  getOri :: o -> Ori
  -- ^ Invert orientation

  invOri :: o -> o

canonicOri :: (Orientable o) => o -> o
-- ^ Transform a Orientable to canonical orientation (LR)
canonicOri o = if getOri o == LR then o else invOri o

unique :: (Hashable a) => [a] -> [a]
-- ^ Eliminate duplicates of a list
unique = HashSet.toList . HashSet.fromList

lPairs :: [a] -> [(a, a)]
-- ^ Get all consecutive pairs of a list
lPairs l = zip l (tail l)

interleavelists :: [a] -> [a] -> [a]
-- ^ Interleave elements of two list
interleavelists l1 l2 = concat . List.transpose $ [l1, l2]

patternError :: String
-- ^ ERROR message to impossible pattern.
patternError = "ERROR: Pattern shouldn't be possible."

traceValue :: (Show a) => a -> a
-- ^ trace for debug
traceValue = traceValueS ""

traceValueS :: (Show a) => String -> a -> a
-- ^ trace for debug, with prefix string
traceValueS str val = trace (str ++ " ---" ++ show val ++ "---") val

evens :: [a] -> [a]
-- ^ take even positions elements of a list (index start in 0)
evens (x : xs) = x : odds xs
evens _ = []

odds :: [a] -> [a]
-- ^ take odd positions elements of a list (index start in 0)
odds (_ : xs) = evens xs
odds _ = []

minWith :: (Foldable t, Ord b) => (a -> b) -> t a -> Maybe a
-- calculate minimum converting values with function
minWith f l =
  if null l
    then Nothing
    else Just $ List.minimumBy (\x y -> compare (f x) (f y)) l

maxWith :: (Foldable t, Ord b) => (a -> b) -> t a -> Maybe a
-- calculate maximum converting values with function
maxWith f l =
  if null l
    then Nothing
    else Just $ List.maximumBy (\x y -> compare (f x) (f y)) l