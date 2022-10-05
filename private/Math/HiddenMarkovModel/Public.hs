{-# LANGUAGE TypeFamilies #-}
module Math.HiddenMarkovModel.Public (
   T(..),
   Discrete, DiscreteTrained,
   Gaussian, GaussianTrained,
   uniform,
   generate,
   generateLabeled,
   probabilitySequence,
   Normalized.logLikelihood,
   Normalized.reveal,

   Trained(..),
   trainSupervised,
   Normalized.trainUnsupervised,
   mergeTrained, finishTraining, trainMany,
   deviation,

   toCSV,
   fromCSV,
   ) where

import qualified Math.HiddenMarkovModel.Public.Distribution as Distr
import qualified Math.HiddenMarkovModel.Normalized as Normalized
import qualified Math.HiddenMarkovModel.CSV as HMMCSV
import Math.HiddenMarkovModel.Private
          (T(..), Trained(..), mergeTrained, toCells, parseCSV)
import Math.HiddenMarkovModel.Utility
          (squareConstant, distance, matrixDistance,
           randomItemProp, normalizeProb, attachOnes)

import qualified Numeric.LAPACK.Matrix.Array as ArrMatrix
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix ((#!))

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as StorableArray
import qualified Data.Array.Comfort.Shape as Shape

import qualified Text.CSV.Lazy.String as CSV

import qualified System.Random as Rnd

import qualified Control.Monad.Exception.Synchronous as ME
import qualified Control.Monad.Trans.State as MS
import qualified Control.Monad.HT as Monad

import qualified Data.NonEmpty as NonEmpty
import Data.Traversable (Traversable, mapAccumL)
import Data.Foldable (Foldable)



type DiscreteTrained symbol sh prob =
         Trained (Distr.Discrete symbol) sh prob
type Discrete symbol sh prob = T (Distr.Discrete symbol) sh prob

type GaussianTrained emiSh stateSh a =
         Trained (Distr.Gaussian emiSh) stateSh a
type Gaussian emiSh stateSh a = T (Distr.Gaussian emiSh) stateSh a


{- |
Create a model with uniform probabilities
for initial vector and transition matrix
given a distribution for the emissions.
You can use this as a starting point for 'Normalized.trainUnsupervised'.
-}
uniform ::
   (Distr.Info typ, Shape.C sh, Class.Real prob) =>
   Distr.T typ sh prob -> T typ sh prob
uniform distr =
   let sh = Distr.statesShape distr
       c = recip $ fromIntegral $ Shape.size sh
   in  Cons {
          initial = Vector.constant sh c,
          transition = squareConstant sh c,
          distribution = distr
       }


probabilitySequence ::
   (Distr.EmissionProb typ, Shape.Indexed sh, Shape.Index sh ~ state,
    Class.Real prob, Distr.Emission typ prob ~ emission, Traversable f) =>
   T typ sh prob -> f (state, emission) -> f prob
probabilitySequence hmm =
   snd
   .
   mapAccumL
      (\index (s, e) ->
         ((transition hmm #!) . flip (,) s,
          index s * Distr.emissionStateProb (distribution hmm) e s))
      (initial hmm StorableArray.!)

generate ::
   (Distr.Generate typ, Shape.Indexed sh, Class.Real prob,
    Rnd.RandomGen g, Rnd.Random prob, Distr.Emission typ prob ~ emission) =>
   T typ sh prob -> g -> [emission]
generate hmm = map snd . generateLabeled hmm

generateLabeled ::
   (Distr.Generate typ, Shape.Indexed sh, Shape.Index sh ~ state,
    Rnd.RandomGen g, Rnd.Random prob,
    Class.Real prob, Distr.Emission typ prob ~ emission) =>
   T typ sh prob -> g -> [(state, emission)]
generateLabeled hmm =
   MS.evalState $
   flip MS.evalStateT (initial hmm) $
   Monad.repeat $ MS.StateT $ \v0 -> do
      s <-
         randomItemProp $
         zip (Shape.indices $ StorableArray.shape v0) (Vector.toList v0)
      x <- Distr.generate (distribution hmm) s
      return ((s, x), Matrix.takeColumn (transition hmm) s)



{- |
Contribute a manually labeled emission sequence to a HMM training.
-}
trainSupervised ::
   (Distr.Estimate typ, Shape.Indexed sh, Shape.Index sh ~ state,
    Class.Real prob, Distr.Emission typ prob ~ emission) =>
   sh -> NonEmpty.T [] (state, emission) -> Trained typ sh prob
trainSupervised sh xs =
   let getState (s, _x) = s
   in  Trained {
          trainedInitial = Vector.unit sh $ getState $ NonEmpty.head xs,
          trainedTransition =
             Matrix.transpose $ ArrMatrix.fromVector $
             StorableArray.accumulate (+)
                (ArrMatrix.toVector $ squareConstant sh 0) $
             attachOnes $ NonEmpty.mapAdjacent (,) $ fmap getState xs,
          trainedDistribution = Distr.accumulateEmissions sh xs
       }

finishTraining ::
   (Distr.Estimate typ, Shape.C sh, Eq sh, Class.Real prob) =>
   Trained typ sh prob -> T typ sh prob
finishTraining hmm =
   Cons {
      initial = normalizeProb $ trainedInitial hmm,
      transition = normalizeProbColumns $ trainedTransition hmm,
      distribution = Distr.normalize $ trainedDistribution hmm
   }

normalizeProbColumns ::
   (Shape.C sh, Eq sh, Class.Real a) => Matrix.Square sh a -> Matrix.Square sh a
normalizeProbColumns m =
   Matrix.scaleColumns (StorableArray.map recip (Matrix.columnSums m)) m

trainMany ::
   (Distr.Estimate typ, Shape.C sh, Eq sh, Class.Real prob, Foldable f) =>
   (trainingData -> Trained typ sh prob) ->
   NonEmpty.T f trainingData -> T typ sh prob
trainMany train = finishTraining . NonEmpty.foldl1Map mergeTrained train





{- |
Compute maximum deviation between initial and transition probabilities.
You can use this as abort criterion for unsupervised training.
We omit computation of differences between the emission probabilities.
This simplifies matters a lot and
should suffice for defining an abort criterion.
-}
deviation ::
   (Shape.C sh, Eq sh, Class.Real prob) =>
   T typ sh prob -> T typ sh prob -> prob
deviation hmm0 hmm1 =
   distance (initial hmm0) (initial hmm1)
   `max`
   matrixDistance (transition hmm0) (transition hmm1)


toCSV ::
   (Distr.ToCSV typ, Shape.Indexed sh, Class.Real prob, Show prob) =>
   T typ sh prob -> String
toCSV hmm =
   CSV.ppCSVTable $ snd $ CSV.toCSVTable $ HMMCSV.padTable "" $ toCells hmm

fromCSV ::
   (Distr.FromCSV typ, Shape.Indexed sh, Eq sh, Class.Real prob, Read prob) =>
   (Int -> sh) -> String -> ME.Exceptional String (T typ sh prob)
fromCSV makeShape =
   MS.evalStateT (parseCSV makeShape) . map HMMCSV.fixShortRow . CSV.parseCSV
