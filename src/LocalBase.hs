{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE ImportQualifiedPost #-}

module LocalBase where

-- Module      : LocalBase
-- Description : Auxiliar functions.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com

import Data.HashSet qualified as HashSet
import Data.Hashable (Hashable)
import Data.List (sortBy)
import Data.List qualified as List
import Debug.Trace (trace)

newtype Dist = Dist Int deriving newtype (Eq, Show, Read)

data Ori = LR | RL deriving (Eq, Show)

class Orientable o where
  -- Get orientation
  getOri :: o -> Ori

  -- Invert orientation
  invOri :: o -> o

  -- Check if there is two orientations
  hasOri :: o -> Bool

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

indexError :: String
-- ^ ERROR message for index out of bounds.
indexError = "ERROR: Index out of bounds."

inputError :: String
-- ^ ERROR message to incorrect input.
inputError = "ERROR: Input is invalid."

logicError :: String
-- ^ ERROR message for error in the logic or implementation.
logicError = "ERROR: There is something wrong with the algorithm or the implementation does not correspond to it."

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

replace :: Int -> (a -> a) -> [a] -> [a]
-- ^ replace element in a list. If you are using this function a lot, maybe you should use a different data structure.
replace i f xs =
  case splitAt i xs of
    (_, []) -> error indexError
    (before, e : after) -> before ++ [f e] ++ after

rotateL :: Int -> [a] -> [a]
-- ^ rotate the given number of positions to the left
rotateL = drop <> take

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

sortWith :: Ord b => (a -> b) -> [a] -> [a]
-- sort list converting values with function
sortWith f = sortBy (\x y -> compare (f x) (f y))

xOr :: Bool -> Bool -> Bool
-- ^ XOR boolean operator
xOr p q = (p || q) && not (p && q)
