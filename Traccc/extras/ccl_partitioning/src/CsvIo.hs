{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE LambdaCase #-}

module CsvIo (getCsvCells) where

import System.Exit (exitFailure)
import System.IO
import Data.ByteString (ByteString, hGetSome, empty)
import Data.Csv.Incremental
import Data (Cell)

feed :: (ByteString -> Parser Cell) -> Handle -> IO (Parser Cell)
feed k csvFile = do
  hIsEOF csvFile >>= \case
    True  -> return $ k empty
    False -> k <$> hGetSome csvFile 4096

getCsvCells :: String -> IO [Either String Cell]
getCsvCells fn = do
  withFile fn ReadMode $ \ csvFile -> do
    let loop !_ (Fail _ errMsg) = do putStrLn errMsg; exitFailure
        loop acc (Many rs k)    = loop (acc ++ rs) =<< feed k csvFile
        loop acc (Done rs)      = do return (acc ++ rs)
    loop [] (decode HasHeader)
