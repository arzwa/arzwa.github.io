# This file was generated, do not modify it. # hide
# Note: in julia one can define short functions using the notation
# `f(x) = <expression in x>`, without using the `function` keyword
Pmatrix(p) = [1-p p/3 p/3 p/3 ;
              p/3 1-p p/3 p/3 ;
              p/3 p/3 1-p p/3 ;
              p/3 p/3 p/3 1-p ]
P = Pmatrix(0.2)