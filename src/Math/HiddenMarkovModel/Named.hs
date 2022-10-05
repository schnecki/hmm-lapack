{-# LANGUAGE TypeFamilies #-}
module Math.HiddenMarkovModel.Named (
   T(..),
   Discrete,
   Gaussian,
   fromModelAndNames,
   toCSV,
   fromCSV,
   ) where

import qualified Math.HiddenMarkovModel.Distribution as Distr
import qualified Math.HiddenMarkovModel.Private as HMM
import qualified Math.HiddenMarkovModel.CSV as HMMCSV
import Math.HiddenMarkovModel.Utility (attachOnes, vectorDim)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as StorableArray
import qualified Data.Array.Comfort.Boxed as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Boxed (Array)

import qualified Text.CSV.Lazy.String as CSV
import Text.Printf (printf)

import qualified Control.Monad.Exception.Synchronous as ME
import qualified Control.Monad.Trans.State as MS
import Control.DeepSeq (NFData, rnf)
import Foreign.Storable (Storable)

import qualified Data.Map as Map
import qualified Data.List as List
import Data.Tuple.HT (swap)
import Data.Map (Map)


{- |
A Hidden Markov Model with names for each state.

Although 'nameFromStateMap' and 'stateFromNameMap' are exported
you must be careful to keep them consistent when you alter them.
-}
data T typ sh ix prob =
   Cons {
      model :: HMM.T typ sh prob,
      nameFromStateMap :: Array sh String,
      stateFromNameMap :: Map String ix
   }
   deriving (Show)

type Simple typ sh prob = T typ sh (Shape.Index sh) prob
type Discrete symbol stateSh prob =
      Simple (Distr.Discrete symbol) stateSh prob
type Gaussian emiSh stateSh a =
      Simple (Distr.Gaussian emiSh) stateSh a


instance
   (Distr.NFData typ, NFData sh, NFData ix, NFData prob,
    Shape.C sh, Storable prob) =>
      NFData (T typ sh ix prob) where
   rnf hmm = rnf (model hmm, nameFromStateMap hmm, stateFromNameMap hmm)


fromModelAndNames ::
   (Shape.Indexed sh) =>
   HMM.T typ sh prob -> [String] -> Simple typ sh prob
fromModelAndNames md names =
   let m = Array.fromList (StorableArray.shape $ HMM.initial md) names
   in  Cons {
          model = md,
          nameFromStateMap = m,
          stateFromNameMap = inverseMap m
       }

inverseMap ::
   (Shape.Indexed sh, Shape.Index sh ~ ix) => Array sh String -> Map String ix
inverseMap =
   Map.fromListWith (error "duplicate label") .
   map swap . Array.toAssociations


toCSV ::
   (Distr.ToCSV typ, Shape.Indexed sh, Class.Real prob, Show prob) =>
   Simple typ sh prob -> String
toCSV hmm =
   CSV.ppCSVTable $ snd $ CSV.toCSVTable $ HMMCSV.padTable "" $
      Array.toList (nameFromStateMap hmm) : HMM.toCells (model hmm)

fromCSV ::
   (Distr.FromCSV typ, Shape.Indexed stateSh, Eq stateSh,
    Class.Real prob, Read prob) =>
   (Int -> stateSh) ->
   String -> ME.Exceptional String (Simple typ stateSh prob)
fromCSV makeShape =
   MS.evalStateT (parseCSV makeShape) . map HMMCSV.fixShortRow . CSV.parseCSV

parseCSV ::
   (Distr.FromCSV typ, Shape.Indexed stateSh, Eq stateSh,
    Class.Real prob, Read prob) =>
   (Int -> stateSh) -> HMMCSV.CSVParser (Simple typ stateSh prob)
parseCSV makeShape = do
   names <- HMMCSV.parseStringList =<< HMMCSV.getRow
   let duplicateNames =
         Map.keys $ Map.filter (> (1::Int)) $
         Map.fromListWith (+) $ attachOnes names
    in HMMCSV.assert (null duplicateNames) $
          "duplicate names: " ++ List.intercalate ", " duplicateNames
   md <- HMM.parseCSV makeShape
   let n = length names
       m = vectorDim (HMM.initial md)
    in HMMCSV.assert (n == m) $
          printf "got %d state names for %d states" n m
   return $ fromModelAndNames md names
