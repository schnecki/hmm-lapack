{-# LANGUAGE TypeFamilies #-}
module Math.HiddenMarkovModel.Example.TrafficLightPrivate where

import qualified Math.HiddenMarkovModel.Public as HMM
import qualified Math.HiddenMarkovModel.Public.Distribution as Distr
import Math.HiddenMarkovModel.Utility (normalizeProb, squareFromLists)

import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Vector (Vector)

import qualified Data.Array.Comfort.Shape as Shape

import Text.Read.HT (maybeRead)

import Control.DeepSeq (NFData(rnf))

import qualified Data.NonEmpty as NonEmpty
import qualified Data.List.HT as ListHT
import Data.NonEmpty ((!:))


{- $setup
>>> import qualified Data.NonEmpty as NonEmpty
>>> import Control.DeepSeq (deepseq)
>>>
>>> verifyRevelations :: HMM -> [Bool]
>>> verifyRevelations hmm_ =
>>>    map (verifyRevelation hmm_) (NonEmpty.flatten labeledSequences)
-}


data Color = Red | Yellow | Green
   deriving (Eq, Ord, Enum, Show, Read)

instance NFData Color where
   rnf Red = ()
   rnf _ = ()

{- |
Using 'show' and 'read' is not always a good choice
since they must format and parse Haskell expressions
which is not of much use to the outside world.
-}
instance Distr.CSVSymbol Color where
   cellFromSymbol = show
   symbolFromCell = maybeRead


data State = StateRed | StateYellowRG | StateGreen | StateYellowGR
   deriving (Eq, Ord, Enum, Bounded)

type StateSet = Shape.Enumeration State

stateSet :: StateSet
stateSet = Shape.Enumeration


type HMM = HMM.Discrete Color StateSet Double

{- |
>>> verifyRevelations hmm
[True,True]
-}
hmm :: HMM
hmm =
   HMM.Cons {
      HMM.initial = normalizeProb $ stateVector 2 1 2 1,
      HMM.transition =
         squareFromLists stateSet $
            stateVector 0.8 0.0 0.0 0.2 :
            stateVector 0.2 0.8 0.0 0.0 :
            stateVector 0.0 0.2 0.8 0.0 :
            stateVector 0.0 0.0 0.2 0.8 :
            [],
      HMM.distribution =
         Distr.discreteFromList $
            (Red,    stateVector 1 0 0 0) !:
            (Yellow, stateVector 0 1 0 1) :
            (Green,  stateVector 0 0 1 0) :
            []
   }


{- |
>>> verifyRevelations hmmDisturbed
[True,True]
-}
hmmDisturbed :: HMM
hmmDisturbed =
   HMM.Cons {
      HMM.initial = normalizeProb $ stateVector 1 1 1 1,
      HMM.transition =
         squareFromLists stateSet $
            stateVector 0.3 0.2 0.2 0.3 :
            stateVector 0.3 0.3 0.2 0.2 :
            stateVector 0.2 0.3 0.3 0.2 :
            stateVector 0.2 0.2 0.3 0.3 :
            [],
      HMM.distribution =
         Distr.discreteFromList $
            (Red,    stateVector 0.6 0.2 0.2 0.2) !:
            (Yellow, stateVector 0.2 0.6 0.2 0.6) :
            (Green,  stateVector 0.2 0.2 0.6 0.2) :
            []
   }

stateVector :: Double -> Double -> Double -> Double -> Vector StateSet Double
stateVector x0 x1 x2 x3 = Vector.fromList stateSet [x0,x1,x2,x3]


red, yellowRG, green, yellowGR :: (State, Color)
red      = (StateRed, Red)
yellowRG = (StateYellowRG, Yellow)
green    = (StateGreen, Green)
yellowGR = (StateYellowGR, Yellow)

labeledSequences :: NonEmpty.T [] (NonEmpty.T [] (State, Color))
labeledSequences =
   (red !: red : red : red :
    yellowRG : yellowRG :
    green : green : green : green : green :
    yellowGR :
    red : red : red :
    []) !:
   (green !: green : green :
    yellowGR :
    red : red : red : red :
    yellowRG :
    green : green : green : green : green :
    yellowGR : yellowGR :
    []) :
   []

{- |
Construct a Hidden Markov model by watching a set
of manually created sequences of emissions and according states.

>>> verifyRevelations hmmTrainedSupervised
[True,True]
-}
hmmTrainedSupervised :: HMM
hmmTrainedSupervised =
   HMM.trainMany (HMM.trainSupervised stateSet) labeledSequences


stateSequences :: NonEmpty.T [] (NonEmpty.T [] Color)
stateSequences = fmap (fmap snd) labeledSequences

{- |
Construct a Hidden Markov model starting from a known model
and a set of sequences that contain only the emissions, but no states.

>>> verifyRevelations hmmTrainedUnsupervised
[True,True]
-}
hmmTrainedUnsupervised :: HMM
hmmTrainedUnsupervised =
   HMM.trainMany (HMM.trainUnsupervised hmm) stateSequences

{- |
Repeat unsupervised training until convergence.

prop> deepseq hmmIterativelyTrained True
-}
hmmIterativelyTrained :: HMM
hmmIterativelyTrained =
   snd $ head $ dropWhile fst $
   ListHT.mapAdjacent (\hmm0 hmm1 -> (HMM.deviation hmm0 hmm1 > 1e-5, hmm1)) $
   iterate
      (flip HMM.trainMany stateSequences . HMM.trainUnsupervised)
      hmmDisturbed


verifyRevelation :: HMM -> NonEmpty.T [] (State, Color) -> Bool
verifyRevelation model xs =
   fmap fst xs == HMM.reveal model (fmap snd xs)
