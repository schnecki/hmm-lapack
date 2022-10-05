{-# LANGUAGE TypeFamilies #-}
{- |
This module provides a simple way to train
the transition matrix and initial probability vector
using simple patterns of state sequences.

You may create a trained model using semigroup combinators like this:

> example :: HMM.DiscreteTrained Char (ShapeStatic.ZeroBased TypeNum.U2) Double
> example =
>    let a = atom FL.i0
>        b = atom FL.i1
>        distr =
>           Distr.DiscreteTrained $ Map.fromList $
>           ('a', ShapeStatic.vector $ 1!:2!:FL.end) :
>           ('b', ShapeStatic.vector $ 4!:3!:FL.end) :
>           ('c', ShapeStatic.vector $ 0!:1!:FL.end) :
>           []
>    in finish (ShapeStatic.ZeroBased Proxy) distr $
>       replicate 5 $ replicate 10 a <> replicate 20 b
-}
module Math.HiddenMarkovModel.Pattern (
   T,
   atom,
   append,
   replicate,
   finish,
   ) where

import qualified Math.HiddenMarkovModel.Distribution as Distr
import qualified Math.HiddenMarkovModel as HMM
import Math.HiddenMarkovModel.Private (Trained(..))
import Math.HiddenMarkovModel.Utility (squareConstant)

import qualified Numeric.LAPACK.Matrix.Array as ArrMatrix
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as StorableArray
import qualified Data.Array.Comfort.Shape.Static as ShapeStatic
import qualified Data.Array.Comfort.Shape as Shape

import qualified Data.FixedLength as FL
import Data.FixedLength ((!:))

import qualified Type.Data.Num.Unary.Literal as TypeNum

import qualified Data.NonEmpty.Map as NonEmptyMap
import qualified Data.NonEmpty as NonEmpty
import Data.Semigroup (Semigroup, (<>), stimes)

import Prelude hiding (replicate)


newtype T sh prob =
   Cons (sh -> (Shape.Index sh, Matrix.Square sh prob, Shape.Index sh))

atom ::
   (Shape.Indexed sh, Shape.Index sh ~ state, Class.Real prob) =>
   state -> T sh prob
atom s = Cons $ \sh -> (s, squareConstant sh 0, s)


instance
   (Shape.Indexed sh, Eq sh, Class.Real prob) =>
      Semigroup (T sh prob) where
   (<>) = append
   stimes k = replicate $ fromIntegral k


infixl 5 `append`

append ::
   (Shape.Indexed sh, Eq sh, Class.Real prob) =>
   T sh prob -> T sh prob -> T sh prob
append (Cons f) (Cons g) =
   Cons $ \n ->
      case (f n, g n) of
         ((sai, ma, sao), (sbi, mb, sbo)) ->
            (sai, increment (sbi,sao) 1 $ Matrix.add ma mb, sbo)

replicate ::
   (Shape.Indexed sh, Class.Real prob) => Int -> T sh prob -> T sh prob
replicate ki (Cons f) =
   Cons $ \sh ->
      case f sh of
         (si, m, so) ->
            let k = fromIntegral ki
            in  (si, increment (si,so) (k-1) $ Matrix.scale k m, so)

increment ::
   (Shape.Indexed sh, Shape.Index sh ~ state, Class.Real a) =>
   (state, state) -> a -> Matrix.Square sh a -> Matrix.Square sh a
increment (i,j) x =
   ArrMatrix.lift1 $ flip (StorableArray.accumulate (+)) [((i,j), x)]


finish ::
   (Distr.Info typ, Shape.Indexed sh, Class.Real prob) =>
   Distr.Trained typ sh prob -> T sh prob -> Trained typ sh prob
finish tdistr (Cons f) =
   let sh = Distr.statesShapeTrained tdistr
       (si, m, _so) = f sh
   in Trained {
         trainedInitial = Vector.unit sh si,
         trainedTransition = m,
         trainedDistribution = tdistr
      }


_example :: HMM.DiscreteTrained Char (ShapeStatic.ZeroBased TypeNum.U2) Double
_example =
   let a = atom FL.i0
       b = atom FL.i1
       distr =
          Distr.DiscreteTrained $ NonEmptyMap.fromList $
          ('a', ShapeStatic.vector $ 1!:2!:FL.end) NonEmpty.!:
          ('b', ShapeStatic.vector $ 4!:3!:FL.end) :
          ('c', ShapeStatic.vector $ 0!:1!:FL.end) :
          []
   in finish distr $ replicate 5 $ replicate 10 a <> replicate 20 b
