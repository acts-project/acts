{-# LANGUAGE DeriveGeneric #-}

module Data (Cell(..), Partition(..)) where

import GHC.Generics
import Data.Csv (FromRecord, ToRecord)

data Cell = Cell
  { moduleId :: Integer
  , measurementId :: Integer
  , channel0 :: Integer
  , channel1 :: Integer
  , timestamp :: Float
  , value :: Float
  } deriving (Show, Eq, Generic)

data Partition = Partition
  { start :: Integer
  , psize :: Integer
  } deriving (Show, Eq)

instance FromRecord Cell
instance ToRecord Cell
