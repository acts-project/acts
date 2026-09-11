import System.Exit (exitFailure)
import Data.Maybe (isNothing, fromJust)
import Data.Either (partitionEithers)
import Control.Monad (forM_)
import System.Environment
import CsvIo (getCsvCells)
import Data (Cell(..), Partition(..))

makePartitions :: [Cell] -> [Partition]
makePartitions [] = []
makePartitions acs = go acs 0 0 0 Nothing []
  where
    go :: [Cell] -> Integer -> Integer -> Integer -> Maybe Cell -> [Partition] -> [Partition]
    go [] _ s l _ p = p ++ [Partition { start = s, psize = l }]
    go (c:cs) i s l lc p
      | isNothing lc = go cs (i + 1) 0 1 (Just c) p
      | moduleId c /= moduleId (fromJust lc) = go cs (i + 1) i 1 (Just c) (p ++ [Partition { start = s, psize = l }])
      | otherwise = go cs (i + 1) s (l + 1) (Just c) p

packPartitionsTrivial :: [Partition] -> [[Partition]]
packPartitionsTrivial [] = []
packPartitionsTrivial a = [a]

packPartitionsClassic :: Integer -> [Partition] -> [[Partition]]
packPartitionsClassic bs aps = go aps 0 [] []
  where
    go :: [Partition] -> Integer -> [Partition] -> [[Partition]] -> [[Partition]]
    go [] _ c r
      | null c = r
      | otherwise = r ++ [c]
    go (p:ps) n c r
      | ((psize p) + n) >= bs = go ps 0 [] (r ++ [(c ++ [p])])
      | otherwise = go ps (n + (psize p)) (c ++ [p]) r

printPackingStats :: [[Partition]] -> IO ()
printPackingStats p = do
  putStrLn ("Total partitions: " ++ (show . length $ concatPartitions))
  putStrLn ("Total cells:      " ++ (show numCells))
  putStrLn ("Total bins:       " ++ (show numBins))
  putStrLn ("Smallest bin:     " ++ (show smallestBin))
  putStrLn ("Largest bin:      " ++ (show largestBin))
  putStrLn ("Inefficiency:     " ++ (show ineffRatio))
  putStrLn ("Bins:             " ++ (show binCells))
  where
    numBins = toInteger (length p)
    concatPartitions = concat p
    numCells = sum . map psize $ concatPartitions
    binCells = map (sum . map psize) $ p
    smallestBin = minimum binCells
    largestBin = maximum binCells
    ineffRatio = ((fromInteger (numBins * largestBin)) / (fromInteger (numCells))) :: Float


main :: IO ()
main = do
  args <- getArgs;
  let fileName = args !! 0
  putStrLn ("Reading cells from " ++ fileName)
  errorCells <- getCsvCells fileName
  let (errors, cells) = partitionEithers errorCells
  if not (null errors) then
    do
      putStrLn "Some cells failed to parse";
      exitFailure
    else
      putStrLn "Successfully read all cells!"
  let partitions = makePartitions cells
  let packingAlgs =
        [ ("Trivial packing", packPartitionsTrivial)
        , ("Classic packing (n = 1024)", (packPartitionsClassic 1024))
        , ("Classic packing (n = 2048)", (packPartitionsClassic 2048))
        ]
  forM_ packingAlgs (\(name, fun) -> do
    putStrLn ("\n" ++ name ++ ":")
    printPackingStats (fun partitions))
