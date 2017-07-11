#lang scribble/manual
@(require scribble/examples)
@(require "pr-math.rkt")
@require[@for-label[octavian
                    racket]]

@setup-math
@(define octavian-eval (make-base-eval))
@(octavian-eval '(require octavian))

@title{octavian: scientific computing methods}
@author{Luis Osa}

@defmodule[octavian]

The package @italic{octavian} provides scientific computing methods resembling those
found in Octave & MATLAB.

@table-of-contents[]

@include-section["nonlinear.scrbl"]
@include-section["approximation.scrbl"]