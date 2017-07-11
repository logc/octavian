#lang scribble/manual
@(require scribble/examples)
@(require "pr-math.rkt")
@require[@for-label[octavian
                    racket]]

@setup-math
@(define octavian-eval (make-base-eval))
@(octavian-eval '(require octavian))

@title{Approximation of functions and data}


@defproc[(polyfit [xs (listof real?)] [ys (listof real?)] [grade real?])
         (col-matrix?)]
fits a function @${f(x)}, described by a list of coordinates @racket[xs] and a
list of values @racket[ys], with a polynomial of up to a specified
@racket[grade].  The polynomial @${\tilde{f}} is a good fit for the function in
the sense that it minimizes the sum of squared errors between the approximated
function values and the actual values @${y} at the given coordinates, when
compared with any other approximation using polynomials.

@$${\forall p_m \in \mathbb{P}_m:
\sum_{i=0}^{n}{[y_i - \tilde{f}(x_i)]^2} \leq \sum_{i=0}^{n}{[y_i - p_m(x_i)]^2}}