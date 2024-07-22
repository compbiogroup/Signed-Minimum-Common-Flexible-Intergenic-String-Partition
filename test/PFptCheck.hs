{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE TemplateHaskell #-}

module PFptCheck (tests) where

import Hedgehog
import Hedgehog.Gen qualified as Gen
import Hedgehog.Range qualified as Range
import PFpt

-- prop_graphCostAndPartitionCostAreTheSame :: Property
-- prop_graphCostAndPartitionCostAreTheSame = undefined

tests :: IO Bool
tests = checkSequential $$(discover)
