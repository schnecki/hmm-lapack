module Math.HiddenMarkovModel.Example.SineWavePrivate where

import qualified Math.HiddenMarkovModel.Public as HMM
import qualified Math.HiddenMarkovModel.Public.Distribution as Distr
import Math.HiddenMarkovModel.Utility (normalizeProb, squareFromLists)

import qualified Numeric.LAPACK.Matrix.Hermitian as Hermitian
import qualified Numeric.LAPACK.Matrix.Layout as Layout
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Vector (Vector, singleton)

import qualified Data.Array.Comfort.Boxed as Array
import qualified Data.Array.Comfort.Shape as Shape

import qualified Data.NonEmpty.Class as NonEmptyC
import qualified Data.NonEmpty as NonEmpty
import Data.Function.HT (nest)
import Data.Tuple.HT (mapSnd)


{- $setup
>>> import qualified Data.NonEmpty as NonEmpty
-}


data State = Rising | High | Falling | Low
   deriving (Eq, Ord, Enum, Bounded, Show)

type StateSet = Shape.Enumeration State

stateSet :: StateSet
stateSet = Shape.Enumeration


type HMM = HMM.Gaussian () StateSet Double

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
         Distr.gaussian $ Array.fromList stateSet $
            (singleton   0 , Hermitian.identity Layout.RowMajor ()) :
            (singleton   1 , Hermitian.identity Layout.RowMajor ()) :
            (singleton   0 , Hermitian.identity Layout.RowMajor ()) :
            (singleton (-1), Hermitian.identity Layout.RowMajor ()) :
            []
   }

stateVector :: Double -> Double -> Double -> Double -> Vector StateSet Double
stateVector x0 x1 x2 x3 = Vector.fromList stateSet [x0,x1,x2,x3]

{- |
>>> take 20 $ map fst $ NonEmpty.flatten sineWaveLabeled
[Rising,Rising,High,High,High,Falling,Falling,Falling,Low,Low,Low,Rising,Rising,Rising,Rising,High,High,High,Falling,Falling]
-}
sineWaveLabeled :: NonEmpty.T [] (State, Double)
sineWaveLabeled =
   NonEmpty.mapTail (take 200) $
   fmap (\x -> (toEnum $ mod (floor (x*2/pi+0.5)) 4, sin x)) $
   NonEmptyC.iterate (0.5+) 0

sineWave :: NonEmpty.T [] Double
sineWave = fmap snd sineWaveLabeled

{- |
>>> take 20 $ NonEmpty.flatten revealed
[Rising,Rising,High,High,High,Falling,Falling,Falling,Low,Low,Low,Low,Rising,Rising,Rising,High,High,High,Falling,Falling]
-}
revealed :: NonEmpty.T [] State
revealed = HMM.reveal hmmTrainedSupervised $ fmap singleton sineWave

hmmTrainedSupervised :: HMM
hmmTrainedSupervised =
   HMM.finishTraining $ HMM.trainSupervised stateSet $
   fmap (mapSnd singleton) sineWaveLabeled

hmmTrainedUnsupervised :: HMM
hmmTrainedUnsupervised =
   HMM.finishTraining $ HMM.trainUnsupervised hmm $ fmap singleton sineWave

hmmIterativelyTrained :: HMM
hmmIterativelyTrained =
   nest 100
      (\model ->
         HMM.finishTraining $ HMM.trainUnsupervised model $
         fmap singleton sineWave)
      hmm
