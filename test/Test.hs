module Test (tests) where

import qualified Math.HiddenMarkovModel as HMM
import qualified Math.HiddenMarkovModel.Normalized as Normalized
import qualified Math.HiddenMarkovModel.Private as Priv
import qualified Math.HiddenMarkovModel.Distribution as Distr
import Math.HiddenMarkovModel.Utility
         (squareFromLists, distance, matrixDistance)

import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Vector (Vector)

import qualified Data.Array.Comfort.Shape.Static as ShapeStatic
import qualified Data.Array.Comfort.Shape as Shape

import qualified Data.FixedLength as FL

import qualified Type.Data.Num.Unary.Literal as TypeNum
import Type.Base.Proxy (Proxy(Proxy))

import qualified Test.QuickCheck as QC
import qualified System.Random as Rnd

import qualified Data.NonEmpty.Class as NonEmptyC
import qualified Data.NonEmpty.Map as NonEmptyMap
import qualified Data.NonEmpty as NonEmpty
import qualified Data.Traversable as Trav
import qualified Data.Foldable as Fold
import qualified Data.Map as Map
import Data.NonEmpty ((!:))


type StateSet = ShapeStatic.ZeroBased TypeNum.U4

hmm :: HMM.Discrete Char StateSet Double
hmm =
   HMM.Cons {
      HMM.initial = stateVector 0.1 0.2 0.3 0.4,
      HMM.transition =
         squareFromLists stateSet $
            stateVector 0.7 0.1 0.0 0.2 :
            stateVector 0.1 0.6 0.1 0.0 :
            stateVector 0.1 0.2 0.7 0.0 :
            stateVector 0.1 0.1 0.2 0.8 :
            [],
      HMM.distribution =
         Distr.discreteFromList $
            ('a', stateVector 1 0 0 0) !:
            ('b', stateVector 0 1 0 1) :
            ('c', stateVector 0 0 1 0) :
            []
   }

stateSet :: StateSet
stateSet = ShapeStatic.ZeroBased Proxy

stateVector :: Double -> Double -> Double -> Double -> Vector StateSet Double
stateVector =
   FL.curry
      (ShapeStatic.vector :: FL.T TypeNum.U4 Double -> Vector StateSet Double)


sequ :: NonEmpty.T [] Char
sequ = NonEmpty.cons 'a' $ take 20 (HMM.generate hmm (Rnd.mkStdGen 42))

possibleStates :: Char -> [FL.Index TypeNum.U4]
possibleStates c =
   map fst $ filter snd $
   zip (Shape.indices stateSet) $
   map
      (\p ->
         case p of
            0 -> False
            1 -> True
            _ -> error "invalid emission probability (must be 0 or 1)") $
   Vector.toList $
   case HMM.distribution hmm of Distr.Discrete m -> Matrix.takeRow m c

{- |
Should all be equal.
-}
sequLikelihood :: ((Double, Double), Double, Double, NonEmpty.T [] Double)
sequLikelihood =
   ((Priv.forward hmm sequ, Priv.backward hmm sequ),
    exp $ Normalized.logLikelihood hmm sequ,
    sum $
       map (NonEmpty.product . HMM.probabilitySequence hmm) $
       Trav.mapM (\c -> map (flip (,) c) $ possibleStates c) sequ,
    NonEmptyC.zipWith Vector.dot
       (Priv.alpha hmm sequ)
       (Priv.beta hmm $ NonEmpty.tail sequ))

{- |
Should all be one.
-}
sequLikelihoodNormalized :: NonEmpty.T [] Double
sequLikelihoodNormalized =
   let (calphas,betas) = Normalized.alphaBeta hmm sequ
   in  NonEmptyC.zipWith Vector.dot (fmap snd calphas) betas


{- |
Lists should be equal, but the first list contains one less element.
-}
zetas ::
   ([Vector StateSet Double],
    NonEmpty.T [] (Vector StateSet Double),
    NonEmpty.T [] (Vector StateSet Double))
zetas =
   let (recipLikelihood, alphas, betas) = Priv.alphaBeta hmm sequ
   in  (Priv.zetaFromXi $
           Priv.xiFromAlphaBeta hmm recipLikelihood sequ alphas betas,
        Priv.zetaFromAlphaBeta recipLikelihood alphas betas,
        uncurry Normalized.zetaFromAlphaBeta $
        Normalized.alphaBeta hmm sequ)


{- |
Quick test of zetas - result should be @(True, very small, very small)@.
-}
zetasDiff :: (Bool, Double, Double)
zetasDiff =
   case zetas of
      (z0,z1,z2) ->
         (length z0 == length (NonEmpty.tail z1) &&
          length z0 == length (NonEmpty.tail z2),
          maximum $ zipWith distance z0 $ NonEmpty.init z1,
          NonEmpty.maximum $ NonEmptyC.zipWith distance z1 z2)

{- |
Lists should be equal
-}
xis :: ([Matrix.Square StateSet Double], [Matrix.Square StateSet Double])
xis =
   let (recipLikelihood, alphas, betas) = Priv.alphaBeta hmm sequ
   in  (Priv.xiFromAlphaBeta hmm recipLikelihood sequ alphas betas,
        uncurry (Normalized.xiFromAlphaBeta hmm sequ) $
        Normalized.alphaBeta hmm sequ)

{- |
Quick test of xis - result should be @(True, very small)@.
-}
xisDiff :: (Bool, Double)
xisDiff =
   case xis of
      (x0,x1) ->
         (length x0 == length x1, maximum $ zipWith matrixDistance x0 x1)


reveal :: Bool
reveal =
   Normalized.reveal hmm sequ == Priv.reveal hmm sequ


trainUnsupervised ::
   (HMM.DiscreteTrained Char StateSet Double,
    HMM.DiscreteTrained Char StateSet Double)
trainUnsupervised =
   (Priv.trainUnsupervised hmm sequ,
    Normalized.trainUnsupervised hmm sequ)

trainUnsupervisedDiff :: (Double, Double, (Bool, Double))
trainUnsupervisedDiff =
   case trainUnsupervised of
      (hmm0,hmm1) ->
         (matrixDistance
             (Priv.trainedTransition hmm0) (Priv.trainedTransition hmm1),
          distance
             (Priv.trainedInitial hmm0) (Priv.trainedInitial hmm1),
          case (Priv.trainedDistribution hmm0, Priv.trainedDistribution hmm1) of
             (Distr.DiscreteTrained m0, Distr.DiscreteTrained m1) ->
                (NonEmptyMap.size m0 == NonEmptyMap.size m1,
                 Fold.maximum $
                 Map.intersectionWith distance
                    (NonEmptyMap.flatten m0) (NonEmptyMap.flatten m1)))


allPair :: (a -> Bool, b -> Bool) -> (a,b) -> Bool
allPair (f,g) (a,b) = f a && g b

allTriple :: (a -> Bool, b -> Bool, c -> Bool) -> (a,b,c) -> Bool
allTriple (f,g,h) (a,b,c) = f a && g b && h c

almostZero :: Double -> Bool
almostZero x  =  x < 1e-10

almostOne :: Double -> Bool
almostOne x  =  almostZero $ abs (x-1)

almostEqual :: Double -> Double -> Bool
almostEqual x y  =  almostZero $ abs (x-y)

tests :: [(String, QC.Property)]
tests =
   ("sequLikelihood",
      QC.property $
      case sequLikelihood of
         (forwardBackward, expLog, sumProb, alphaBetas) ->
            allPair (almostEqual sumProb, almostEqual sumProb) forwardBackward
            &&
            almostEqual sumProb expLog
            &&
            length (NonEmpty.tail sequ) == length (NonEmpty.tail alphaBetas)
            &&
            Fold.all (almostEqual sumProb) alphaBetas) :
   ("sequLikelihoodNormalized",
      QC.property $
      length (NonEmpty.tail sequ) ==
         length (NonEmpty.tail sequLikelihoodNormalized)
      &&
      Fold.all almostOne sequLikelihoodNormalized) :
   ("zetasDiff",
      QC.property $ allTriple (id, almostZero, almostZero) zetasDiff) :
   ("xisDiff", QC.property $ allPair (id, almostZero) xisDiff) :
   ("reveal", QC.property reveal) :
   ("trainUnsupervisedDiff",
      QC.property $
      allTriple (almostZero, almostZero, allPair (id, almostZero)) $
      trainUnsupervisedDiff) :
   []
