module Math.HiddenMarkovModel.Example.CirclePrivate where

import qualified Math.HiddenMarkovModel.Public as HMM
import qualified Math.HiddenMarkovModel.Public.Distribution as Distr
import Math.HiddenMarkovModel.Utility
         (normalizeProb, squareFromLists, hermitianFromList)

import qualified Numeric.LAPACK.Matrix.HermitianPositiveDefinite as HermitianPD
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Vector (Vector)

import qualified Data.Array.Comfort.Boxed as Array
import qualified Data.Array.Comfort.Shape as Shape

import qualified System.Random as Rnd

import qualified Control.Monad.Trans.State as MS
import Control.Monad (liftM2, replicateM)

import qualified Data.NonEmpty.Class as NonEmptyC
import qualified Data.NonEmpty as NonEmpty
import Data.Function.HT (nest)
import Data.NonEmpty ((!:))
import Data.Maybe (fromMaybe)


{- $setup
>>> import qualified Math.HiddenMarkovModel as HMM
>>> import qualified Data.NonEmpty as NonEmpty
>>> import Data.Eq.HT (equating)
>>>
>>> checkTraining :: (Int, HMM) -> Bool
>>> checkTraining (maxDiff,hmm_) =
>>>    maxDiff >=
>>>    (length $ filter id $ NonEmpty.flatten $
>>>     NonEmpty.zipWith (/=)
>>>       (HMM.reveal hmm_ circle) (fmap fst circleLabeled))
-}


data State = Q1 | Q2 | Q3 | Q4
   deriving (Eq, Ord, Enum, Bounded, Show)

type StateSet = Shape.Enumeration State

stateSet :: StateSet
stateSet = Shape.Enumeration


data Coordinate = X | Y
   deriving (Eq, Ord, Enum, Bounded)

type CoordinateSet = Shape.Enumeration Coordinate

coordinateSet :: CoordinateSet
coordinateSet = Shape.Enumeration

type HMM = HMM.Gaussian CoordinateSet StateSet Double

{- |
prop> checkTraining (0, hmm)
-}
hmm :: HMM
hmm =
   HMM.Cons {
      HMM.initial = normalizeProb $ Vector.one stateSet,
      HMM.transition =
         squareFromLists stateSet $
            stateVector 0.9 0.0 0.0 0.1 :
            stateVector 0.1 0.9 0.0 0.0 :
            stateVector 0.0 0.1 0.9 0.0 :
            stateVector 0.0 0.0 0.1 0.9 :
            [],
      HMM.distribution =
         let hermitianPD =
               HermitianPD.assurePositiveDefiniteness .
               hermitianFromList coordinateSet
             cov0 = hermitianPD [0.10, -0.09, 0.10]
             cov1 = hermitianPD [0.10,  0.09, 0.10]
         in  Distr.gaussian $ Array.fromList stateSet $
                (Vector.fromList coordinateSet [ 0.5,  0.5], cov0) :
                (Vector.fromList coordinateSet [-0.5,  0.5], cov1) :
                (Vector.fromList coordinateSet [-0.5, -0.5], cov0) :
                (Vector.fromList coordinateSet [ 0.5, -0.5], cov1) :
                []
   }

stateVector :: Double -> Double -> Double -> Double -> Vector StateSet Double
stateVector x0 x1 x2 x3 = Vector.fromList stateSet [x0,x1,x2,x3]

circleLabeled :: NonEmpty.T [] (State, Vector CoordinateSet Double)
circleLabeled =
   NonEmpty.mapTail (take 200) $
   fmap
      (\x ->
         (toEnum $ mod (floor (x*2/pi)) 4,
          Vector.fromList coordinateSet [cos x, sin x])) $
   NonEmptyC.iterate (0.5+) 0

circle :: NonEmpty.T [] (Vector CoordinateSet Double)
circle = fmap snd circleLabeled

{- |
>>> take 20 $ NonEmpty.flatten revealed
[Q1,Q1,Q1,Q1,Q2,Q2,Q2,Q3,Q3,Q3,Q4,Q4,Q4,Q1,Q1,Q1,Q2,Q2,Q2,Q3]

prop> equating (take 1000 . NonEmpty.flatten) revealed $ fmap fst circleLabeled
-}
revealed :: NonEmpty.T [] State
revealed = HMM.reveal hmm circle

{- |
Sample multivariate normal distribution and reconstruct it from the samples.
You should obtain the same parameters.
-}
reconstructDistribution :: HMM.Gaussian CoordinateSet () Double
reconstructDistribution =
   let gen = Distr.generate (HMM.distribution hmm) Q1
   in  HMM.finishTraining $ HMM.trainSupervised () $ fmap ((,) ()) $
       flip MS.evalState (Rnd.mkStdGen 23) $
       liftM2 (!:) gen $ replicateM 1000 gen

{- |
Generate labeled emission sequences
and use them for supervised training.

prop> checkTraining (0, reconstructModel)
-}
reconstructModel :: HMM
reconstructModel =
   HMM.trainMany (HMM.trainSupervised stateSet) $
   fmap
      (\seed ->
         fromMaybe (error "empty generated sequence") $ NonEmpty.fetch $
         take 1000 $ HMM.generateLabeled hmm $ Rnd.mkStdGen seed)
      (23 !: take 42 [24..])


{- |
prop> checkTraining (0, hmmTrainedSupervised)
-}
hmmTrainedSupervised :: HMM
hmmTrainedSupervised =
   HMM.finishTraining $ HMM.trainSupervised stateSet circleLabeled

{- |
prop> checkTraining (0, hmmTrainedUnsupervised)
-}
hmmTrainedUnsupervised :: HMM
hmmTrainedUnsupervised =
   HMM.finishTraining $ HMM.trainUnsupervised hmm circle

{- |
prop> checkTraining (40, hmmIterativelyTrained)
-}
hmmIterativelyTrained :: HMM
hmmIterativelyTrained =
   nest 100
      (HMM.finishTraining . flip HMM.trainUnsupervised circle)
      hmm
