import GenomesCheck as GC
import PartitionCheck as PC

main :: IO ()
main = do
  ans <- fmap and . sequence $
    [ return True
    , GC.tests
    , PC.tests
    ]
  print ans
