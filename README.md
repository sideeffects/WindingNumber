# Fast Winding Numbers for Soups

Implementation of the _ACM SIGGRAPH_ 2018 paper, 

"Fast Winding Numbers for Soups and Clouds" 

Gavin Barill¹, Neil Dickson², Ryan Schmidt³, David I.W. Levin¹, Alec Jacobson¹

¹University of Toronto, ²SideFX, ³Gradient Space


_Note: this implementation is for triangle soups only, not point clouds._

This code, as written, depends on Intel's Threading Building Blocks (TBB) library for parallelism, but it should be fairly easy to change it to use any other means of threading, since it only uses parallel for loops with simple partitioning.
