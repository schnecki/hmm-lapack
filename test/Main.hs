module Main where

import Test (tests)

import qualified Test.DocTest.Driver as DocTest
import qualified Test.Main as TestMain


main :: IO ()
main = DocTest.run $ do
   mapM_
      (\(name,prop) ->
         DocTest.printPrefix (name ++ ": ") >> DocTest.property prop)
      tests
   TestMain.main
