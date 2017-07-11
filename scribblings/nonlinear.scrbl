#lang scribble/manual
@(require scribble/examples)
@(require "pr-math.rkt")
@require[@for-label[octavian
                    racket]]

@setup-math
@(define octavian-eval (make-base-eval))
@(octavian-eval '(require octavian))

@title{Non-linear equations}

@defmodule[octavian/nonlinear]

Procedures defined for non-linear equations are algorithms that find
@italic{zeroes} or @italic{fixed points} of a function. Functions are expected
to be real, i.e. they map real variables to real values. Zeroes of a function
are those values @${x} where: @${f(x) = 0}. Fixed points of a
function are values @${x} where: @${f(x) = x}.

@defproc[(bisection [f (-> real? real?)] [a real?] [b real?] [𝛆 real?]
          [nmax (k_min a b 𝛆)])
         (real-in a b)]
finds a zero of function @${f} in interval @${[a, b]}, with tolerance
@${𝛆} in the result.  Optionally, pass a maximum number of iterations
@italic{nmax}.  There is a procedure, @racket[k_min], to obtain a number of steps
@italic{k_min} that are enough to reach the zero; by default, this formula is
used if no @italic{nmax} is supplied.

@examples[#:eval octavian-eval
(let ([𝛆 1e-5]
      [f (lambda (x) x)])
    (bisection f 0 1 𝛆))
]


@defproc[(k_min [a real?] [b real?] [𝛆 real?])
          (integer?)]
computes the minimum number of steps that are necessary to find a zero in the
interval @${[a, b]} with tolerance @${𝛆} by the procedure @racket[bisection].


@defproc[(newton [f (-> real? real?)] [df (-> real? real?)] [x_0 real?] [𝛆 real?]
          [nmax integer?])
         (real?)]
fnd zeroes of function @${f} using Newton's method, which requires the
derivative function @${df} and a starting point @${x_0}.  The
function also takes a tolerance @${𝛆} and a maximum number of iterations
@italic{nmax}.


@defproc[(aitken [𝜙 (-> real? real?)] [x real?] [𝛆 real?] [nmax integer?]
          [niter integer?] [diff real?])
         (real?)]
finds approximations of fixed point @${𝜶} of function @${𝜙} starting
from initial point @${x_0} using Aitken's extrapolation method. The method
stops after a given @italic{nmax} iterations (default: 100), or after the
absolute value of the difference between two consecutive iterations is less
than a given tolerance @${𝛆} (default: 1 * 10e-4).


@defproc[(hoerner [a (listof integer?)] [z real?])
         (real?)]
evaluates efficiently polynomial @italic{a} at point @italic{z}. The polynomial
is expected to be a list of coefficients @${(a_0 ... a_n)}, ordered from lowest
to highest degree.
