#lang scribble/manual
@(require scribble/examples)
@(require "pr-math.rkt")
@require[@for-label[octavian
                    racket]]

@setup-math
@(define octavian-eval (make-base-eval))
@(octavian-eval '(require octavian))

@title{Approximation of functions and data}


@defproc[
(polyfit [xs (listof real?)] [ys (listof real?)] [grade real?])
col-matrix?]{

fits a function @${f(x)}, described by a list of coordinates @racket[xs] and a
list of values @racket[ys], with a polynomial of up to a specified
@racket[grade].  The polynomial @${\widetilde{f}} is a good fit for the function in
the sense that it minimizes the sum of squared errors between the approximated
function values and the actual values @${y} at the given coordinates, when
compared with any other approximation using polynomials.

@$${\forall p_m \in \mathbb{P}_m:
\sum_{i=0}^{n}{[y_i - \widetilde{f}(x_i)]^2} \leq \sum_{i=0}^{n}{[y_i - p_m(x_i)]^2}}

}


@defproc[(chebyshev-nodes [a real?] [b real?] [n real?])
         (listof real?)]{
Returns a list of nodes that are suitable to avoid Runge’s phenomenon when used
to interpolate a function. In particular, in an arbitrary interval [a, b], the so
called Chebyshev-Gauss-Lobatto nodes are the abscissas of equispaced nodes on the
unit semi-circumference. They satisfy the following relationship:

@$${ x_i = \frac{a + b }{2} + \frac{b - a}{2} \hat{x_i} }

where

@$${ \hat{x_i} = -cos(\pi i/n), i = 0,..., n }

}


@defproc[
(barycentric [xs (listof real?)] [ys (listof real?)] [x_1 real?])
real?
]{
computes the value at the abscissae @${x_1} of the polynomial interpolating data
(@${xs}, @${ys}), by using the barycentric formula:

@$${ \Pi_n(x) = \frac{
  \sum_{k=0}^{n}{\frac{w_k}{x - x_k} y_k}
}{
  \sum_{k=0}^{n}{\frac{w_k}{x - x_k}}
}
}

where

@$${w_k = \left(\prod_{\substack{j=0 \\ j \neq k}}^{n}{x_k - x_j}\right)^{-1} , k=0,...,n}
}