{-# LANGUAGE TypeFamilies #-}
module Math.HiddenMarkovModel.Private where

import qualified Math.HiddenMarkovModel.Public.Distribution as Distr
import qualified Math.HiddenMarkovModel.CSV as HMMCSV
import Math.HiddenMarkovModel.Utility (diagonal)

import qualified Numeric.LAPACK.Matrix.Array as ArrMatrix
import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import qualified Numeric.LAPACK.Format as Format
import Numeric.LAPACK.Matrix ((-*#), (##*#), (#*##), (#*|))
import Numeric.LAPACK.Vector (Vector)

import qualified Numeric.Netlib.Class as Class

import Control.DeepSeq (NFData, rnf)
import Control.Applicative ((<$>))

import Foreign.Storable (Storable)

import qualified Data.Array.Comfort.Storable as StorableArray
import qualified Data.Array.Comfort.Shape as Shape

import qualified Data.NonEmpty.Class as NonEmptyC
import qualified Data.NonEmpty as NonEmpty
import qualified Data.Semigroup as Sg
import qualified Data.List as List
import Data.Semigroup ((<>))
import Data.Traversable (Traversable, mapAccumL)
import Data.Tuple.HT (mapFst, mapSnd, swap)


{- |
A Hidden Markov model consists of a number of (hidden) states
and a set of emissions.
There is a vector for the initial probability of each state
and a matrix containing the probability for switching
from one state to another one.
The 'distribution' field points to probability distributions
that associate every state with emissions of different probability.
Famous distribution instances are discrete and Gaussian distributions.
See "Math.HiddenMarkovModel.Distribution" for details.

The transition matrix is transposed
with respect to popular HMM descriptions.
But I think this is the natural orientation, because this way
you can write \"transition matrix times probability column vector\".
-}
data T typ sh prob =
   Cons {
      initial :: Vector sh prob,
      transition :: Matrix.Square sh prob,
      distribution :: Distr.T typ sh prob
   }
   deriving (Show)

instance
   (Distr.NFData typ, NFData sh, Shape.C sh, NFData prob, Storable prob) =>
      NFData (T typ sh prob) where
   rnf (Cons initial_ transition_ distribution_) =
      rnf (initial_, transition_, distribution_)

instance
   (Distr.Format typ, Format.FormatArray sh, Class.Real prob) =>
      Format.Format (T typ sh prob) where
   format fmt (Cons initial_ transition_ distribution_) =
      Format.format fmt (initial_, transition_, distribution_)

mapStatesShape ::
   (Distr.EmissionProb typ, Shape.C sh0, Shape.C sh1) =>
   (sh0 -> sh1) -> T typ sh0 prob -> T typ sh1 prob
mapStatesShape f hmm =
   Cons {
      initial = StorableArray.mapShape f $ initial hmm,
      transition = Square.mapSize f $ transition hmm,
      distribution = Distr.mapStatesShape f $ distribution hmm
   }


emission ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh, Class.Real prob) =>
   T typ sh prob -> Distr.Emission typ prob -> Vector sh prob
emission  =  Distr.emissionProb . distribution


forward ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh, Class.Real prob,
    Distr.Emission typ prob ~ emission, Traversable f) =>
   T typ sh prob -> NonEmpty.T f emission -> prob
forward hmm = Vector.sum . NonEmpty.last . alpha hmm

alpha ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh, Class.Real prob,
    Distr.Emission typ prob ~ emission, Traversable f) =>
   T typ sh prob -> NonEmpty.T f emission -> NonEmpty.T f (Vector sh prob)
alpha hmm (NonEmpty.Cons x xs) =
   NonEmpty.scanl
      (\alphai xi -> Vector.mul (emission hmm xi) (transition hmm #*| alphai))
      (Vector.mul (emission hmm x) (initial hmm))
      xs


backward ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh, Class.Real prob,
    Distr.Emission typ prob ~ emission, Traversable f) =>
   T typ sh prob -> NonEmpty.T f emission -> prob
backward hmm (NonEmpty.Cons x xs) =
   Vector.dot (initial hmm) $
   Vector.mul (emission hmm x) $
   NonEmpty.head $ beta hmm xs

beta ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh, Class.Real prob,
    Distr.Emission typ prob ~ emission, Traversable f) =>
   T typ sh prob -> f emission -> NonEmpty.T f (Vector sh prob)
beta hmm =
   NonEmpty.scanr
      (\xi betai -> Vector.mul (emission hmm xi) betai -*# transition hmm)
      (Vector.one $ StorableArray.shape $ initial hmm)


alphaBeta ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh, Class.Real prob,
    Distr.Emission typ prob ~ emission, Traversable f) =>
   T typ sh prob ->
   NonEmpty.T f emission ->
   (prob, NonEmpty.T f (Vector sh prob), NonEmpty.T f (Vector sh prob))
alphaBeta hmm xs =
   let alphas = alpha hmm xs
       betas = beta hmm $ NonEmpty.tail xs
       recipLikelihood = recip $ Vector.sum $ NonEmpty.last alphas
   in  (recipLikelihood, alphas, betas)



biscaleTransition ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh, Class.Real prob) =>
   T typ sh prob -> Distr.Emission typ prob ->
   Vector sh prob -> Vector sh prob -> Matrix.Square sh prob
biscaleTransition hmm x alpha0 beta1 =
   (diagonal (Vector.mul (emission hmm x) beta1)
    #*##
    transition hmm)
   ##*#
   diagonal alpha0

xiFromAlphaBeta ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh, Class.Real prob,
    Distr.Emission typ prob ~ emission) =>
   T typ sh prob -> prob ->
   NonEmpty.T [] emission ->
   NonEmpty.T [] (Vector sh prob) ->
   NonEmpty.T [] (Vector sh prob) ->
   [Matrix.Square sh prob]
xiFromAlphaBeta hmm recipLikelihood xs alphas betas =
   zipWith3
      (\x alpha0 beta1 ->
         Matrix.scale recipLikelihood $
         biscaleTransition hmm x alpha0 beta1)
      (NonEmpty.tail xs)
      (NonEmpty.init alphas)
      (NonEmpty.tail betas)

zetaFromXi ::
   (Shape.C sh, Eq sh, Class.Real prob) =>
   [Matrix.Square sh prob] -> [Vector sh prob]
zetaFromXi = map Matrix.columnSums

zetaFromAlphaBeta ::
   (Shape.C sh, Eq sh, Class.Real prob) =>
   prob ->
   NonEmpty.T [] (Vector sh prob) ->
   NonEmpty.T [] (Vector sh prob) ->
   NonEmpty.T [] (Vector sh prob)
zetaFromAlphaBeta recipLikelihood alphas betas =
   fmap (Vector.scale recipLikelihood) $
   NonEmptyC.zipWith Vector.mul alphas betas


{- |
In constrast to Math.HiddenMarkovModel.reveal
this does not normalize the vector.
This is slightly simpler but for long sequences
the product of probabilities might be smaller
than the smallest representable number.
-}
reveal ::
   (Distr.EmissionProb typ, Shape.InvIndexed sh, Eq sh, Shape.Index sh ~ state,
    Distr.Emission typ prob ~ emission, Class.Real prob, Traversable f) =>
   T typ sh prob -> NonEmpty.T f emission -> NonEmpty.T f state
reveal = revealGen id

revealGen ::
   (Distr.EmissionProb typ, Shape.InvIndexed sh, Eq sh, Shape.Index sh ~ state,
    Distr.Emission typ prob ~ emission, Class.Real prob, Traversable f) =>
   (Vector (Shape.Deferred sh) prob -> Vector (Shape.Deferred sh) prob) ->
   T typ sh prob -> NonEmpty.T f emission -> NonEmpty.T f state
revealGen normalize hmm =
   fmap (Shape.revealIndex (StorableArray.shape $ initial hmm)) .
   revealStorable normalize (mapStatesShape Shape.Deferred hmm)

revealStorable ::
   (Distr.EmissionProb typ, Shape.InvIndexed sh, Eq sh,
    Shape.Index sh ~ state, Storable state,
    Distr.Emission typ prob ~ emission, Class.Real prob, Traversable f) =>
   (Vector sh prob -> Vector sh prob) ->
   T typ sh prob -> NonEmpty.T f emission -> NonEmpty.T f state
revealStorable normalize hmm (NonEmpty.Cons x xs) =
   uncurry (NonEmpty.scanr (StorableArray.!)) $
   mapFst (fst . Vector.argAbsMaximum) $
   mapAccumL
      (\alphai xi ->
         swap $ mapSnd (Vector.mul (emission hmm xi)) $
         matrixMaxMul (transition hmm) $ normalize alphai)
      (Vector.mul (emission hmm x) (initial hmm)) xs

matrixMaxMul ::
   (Shape.InvIndexed sh, Eq sh, Shape.Index sh ~ ix, Storable ix,
    Class.Real a) =>
   Matrix.Square sh a -> Vector sh a ->
   (Vector sh ix, Vector sh a)
matrixMaxMul m v = Matrix.rowArgAbsMaximums $ Matrix.scaleColumns v m



{- |
A trained model is a temporary form of a Hidden Markov model
that we need during the training on multiple training sequences.
It allows to collect knowledge over many sequences with 'mergeTrained',
even with mixed supervised and unsupervised training.
You finish the training by converting the trained model
back to a plain modul using 'finishTraining'.

You can create a trained model in three ways:

* supervised training using an emission sequence with associated states,

* unsupervised training using an emission sequence and an existing Hidden Markov Model,

* derive it from state sequence patterns, cf. "Math.HiddenMarkovModel.Pattern".
-}
data Trained typ sh prob =
   Trained {
      trainedInitial :: Vector sh prob,
      trainedTransition :: Matrix.Square sh prob,
      trainedDistribution :: Distr.Trained typ sh prob
   }
   deriving (Show)

instance
   (Distr.NFData typ, NFData sh, Shape.C sh, NFData prob, Storable prob) =>
      NFData (Trained typ sh prob) where
   rnf hmm =
      rnf (trainedInitial hmm, trainedTransition hmm, trainedDistribution hmm)


sumTransitions ::
   (Shape.C sh, Eq sh, Class.Real e) =>
   T typ sh e -> [Matrix.Square sh e] -> Matrix.Square sh e
sumTransitions hmm =
   List.foldl' Matrix.add $
   Matrix.zero $ ArrMatrix.shape $ transition hmm

{- |
Baum-Welch algorithm
-}
trainUnsupervised ::
   (Distr.Estimate typ, Shape.C sh, Eq sh, Class.Real prob,
    Distr.Emission typ prob ~ emission) =>
   T typ sh prob -> NonEmpty.T [] emission -> Trained typ sh prob
trainUnsupervised hmm xs =
   let (recipLikelihood, alphas, betas) = alphaBeta hmm xs
       zetas = zetaFromAlphaBeta recipLikelihood alphas betas
       zeta0 = NonEmpty.head zetas

   in  Trained {
          trainedInitial = zeta0,
          trainedTransition =
             sumTransitions hmm $
             xiFromAlphaBeta hmm recipLikelihood xs alphas betas,
          trainedDistribution =
             Distr.accumulateEmissionVectors $ NonEmptyC.zip xs zetas
       }


mergeTrained ::
   (Distr.Estimate typ, Shape.C sh, Eq sh, Class.Real prob) =>
   Trained typ sh prob -> Trained typ sh prob -> Trained typ sh prob
mergeTrained hmm0 hmm1 =
   Trained {
      trainedInitial = Vector.add (trainedInitial hmm0) (trainedInitial hmm1),
      trainedTransition =
         Matrix.add (trainedTransition hmm0) (trainedTransition hmm1),
      trainedDistribution =
         trainedDistribution hmm0 <> trainedDistribution hmm1
   }

instance
   (Distr.Estimate typ, Shape.C sh, Eq sh, Class.Real prob) =>
      Sg.Semigroup (Trained typ sh prob) where
   (<>) = mergeTrained


toCells ::
   (Distr.ToCSV typ, Shape.Indexed sh, Class.Real prob, Show prob) =>
   T typ sh prob -> [[String]]
toCells hmm =
   (HMMCSV.cellsFromVector $ initial hmm) :
   (HMMCSV.cellsFromSquare $ transition hmm) ++
   [] :
   (Distr.toCells $ distribution hmm)

parseCSV ::
   (Distr.FromCSV typ, Shape.C stateSh, Eq stateSh,
    Class.Real prob, Read prob) =>
   (Int -> stateSh) -> HMMCSV.CSVParser (T typ stateSh prob)
parseCSV makeShape = do
   v <-
      StorableArray.mapShape (makeShape . Shape.zeroBasedSize) <$>
      HMMCSV.parseNonEmptyVectorCells
   let sh = StorableArray.shape v
   m <- HMMCSV.parseSquareMatrixCells sh
   HMMCSV.skipEmptyRow
   distr <- Distr.parseCells sh
   return $ Cons {
      initial = v,
      transition = m,
      distribution = distr
   }
