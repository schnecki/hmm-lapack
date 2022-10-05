module Math.HiddenMarkovModel.CSV where

import Math.HiddenMarkovModel.Utility (vectorDim)

import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix (ShapeInt)
import Numeric.LAPACK.Vector (Vector)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Shape as Shape

import qualified Text.CSV.Lazy.String as CSV
import Text.Read.HT (maybeRead)
import Text.Printf (printf)

import qualified Control.Monad.Exception.Synchronous as ME
import qualified Control.Monad.Trans.Class as MT
import qualified Control.Monad.Trans.State as MS
import Control.Monad.Exception.Synchronous (Exceptional)
import Control.Monad (liftM2, replicateM, unless)

import qualified Data.List.Reverse.StrictElement as Rev
import qualified Data.List.HT as ListHT


cellsFromVector ::
   (Shape.C sh, Show a, Class.Real a) => Vector sh a -> [String]
cellsFromVector = map show . Vector.toList

cellsFromSquare ::
   (Shape.Indexed sh, Show a, Class.Real a) => Matrix.Square sh a -> [[String]]
cellsFromSquare = map (map show . Vector.toList) . Matrix.toRows

padTable :: a -> [[a]] -> [[a]]
padTable x xs =
   let width = maximum (map length xs)
   in  map (ListHT.padRight x width) xs


type CSVParser = MS.StateT CSV.CSVResult (Exceptional String)

assert :: Bool -> String -> CSVParser ()
assert cond msg =
   unless cond $ MT.lift $ ME.throw msg

retrieveShortRow :: CSV.CSVError -> Maybe CSV.CSVRow
retrieveShortRow err =
   case err of
      CSV.IncorrectRow {CSV.csvFields = row} -> Just row
      _ -> Nothing

fixShortRow ::
   Either [CSV.CSVError] CSV.CSVRow -> Either [CSV.CSVError] CSV.CSVRow
fixShortRow erow =
   case erow of
      Left errs ->
         case ListHT.partitionMaybe retrieveShortRow errs of
            ([row], []) -> Right row
            _ -> Left errs
      _ -> erow

maybeGetRow :: CSVParser (Maybe CSV.CSVRow)
maybeGetRow = do
   csv0 <- MS.get
   case csv0 of
      [] -> return Nothing
      item : csv1 -> do
         MS.put csv1
         case item of
            Right row -> return (Just row)
            Left errors ->
               MT.lift $ ME.throw $ unlines $ map CSV.ppCSVError errors

getRow :: CSVParser CSV.CSVRow
getRow =
   MT.lift . ME.fromMaybe "unexpected end of file" =<< maybeGetRow

checkEmptyRow :: CSV.CSVRow -> Exceptional String ()
checkEmptyRow row =
   case filter (not . null . CSV.csvFieldContent) row of
      [] -> return ()
      cell:_ -> ME.throw $ printf "%d: expected empty row" (CSV.csvRowNum cell)

skipEmptyRow :: CSVParser ()
skipEmptyRow  =  MT.lift . checkEmptyRow =<< getRow

manySepUntilEnd :: CSVParser a -> CSVParser [a]
manySepUntilEnd p =
   let go = liftM2 (:) p $ do
          mrow <- maybeGetRow
          case mrow of
             Nothing -> return []
             Just row -> do
                MT.lift $ checkEmptyRow row
                go
   in  go

manyRowsUntilEnd :: (CSV.CSVRow -> CSVParser a) -> CSVParser [a]
manyRowsUntilEnd p =
   let go = do
          mrow <- maybeGetRow
          case mrow of
             Nothing -> return []
             Just row -> liftM2 (:) (p row) go
   in  go

parseVectorCells ::
   (Read a, Class.Real a) =>
   CSVParser (Vector ShapeInt a)
parseVectorCells =
   parseVectorFields =<< getRow

parseVectorFields ::
   (Read a, Class.Real a) =>
   CSV.CSVRow -> CSVParser (Vector ShapeInt a)
parseVectorFields =
   MT.lift . fmap Vector.autoFromList . mapM parseNumberCell .
   Rev.dropWhile (null . CSV.csvFieldContent)

parseNonEmptyVectorCells ::
   (Read a, Class.Real a) =>
   CSVParser (Vector ShapeInt a)
parseNonEmptyVectorCells = do
   v <- parseVectorCells
   assert (vectorDim v > 0) "no data for vector"
   return v

cellContent :: CSV.CSVField -> Exceptional String String
cellContent field =
   case field of
      CSV.CSVFieldError {} -> ME.throw $ CSV.ppCSVField field
      CSV.CSVField { CSV.csvFieldContent = str } -> return str

parseNumberCell :: (Read a) => CSV.CSVField -> Exceptional String a
parseNumberCell field = do
   str <- cellContent field
   ME.fromMaybe (printf "field content \"%s\" is not a number" str) $
      maybeRead str

parseSquareMatrixCells ::
   (Shape.C sh, Read a, Class.Real a) =>
   sh -> CSVParser (Matrix.Square sh a)
parseSquareMatrixCells sh = do
   let n = Shape.size sh
   rows <- replicateM n parseVectorCells
   assert (not $ null rows) "no rows"
   assert (all ((n==) . vectorDim) rows) "inconsistent matrix dimensions"
   return $
      Matrix.reshape (MatrixShape.square MatrixShape.RowMajor sh) $
      Matrix.fromRows (Shape.ZeroBased n) rows

parseStringList :: CSV.CSVRow -> CSVParser [String]
parseStringList =
   MT.lift . mapM cellContent .
   Rev.dropWhile (null . CSV.csvFieldContent)
