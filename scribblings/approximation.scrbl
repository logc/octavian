#lang scribble/manual
@(require scribble/examples)
@(require "pr-math.rkt")
@require[@for-label[octavian
                    racket
                    math/array]]

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


@defproc[
(dft [X (Array A)])
(Array A)
]{
returns the discrete Fourier transform of the array @${X}. Unlike @racket[array-fft], it is
not restricted to arrays of size @${2^n}. For small arrays, a simple naive algorithm is used.

@$${ X_k = \sum_{n=0}^{N-1}{x_n e^{\frac{-2 \pi i k n}{N} }} }

For bigger sizes, if they are a power of 2, then the efficient radix-2 DIT algorithm from
@racket[array-fft] is used. If the bigger algorithm is not a power of 2, then a Cooley-Turkey
Fast Fourier transform is selected. This algorithm does not work for sizes that are a prime
number; for those arrays, an FFT is computed using Rader's algorithm.
}

@defproc[
(idft [X (Array A)])
(Array A)
]{
returns the inverse discrete Fourier transform of the array @${X}. The same algorithms are
selected as for @racket[dft], each one of them in their inverse form.
}


@defproc[
(cooley-turkey-dft [X (Array A)])
(Array A)
]{
computes the discrete FFT of an array @${X} using the Cooley-Turkey algorithm. This is
provided for users who want to select the underlying algorithm of @racket[dft].
}


@defproc[
(rader-dft [X (Array A)])
(Array A)
]{
computes the discrete FFT of an array @${X} using the Rader algorithm. This is provided
for users who want to select the underlying algorithm of @racket[dft].
}


@defproc[
(fftshift [X (Array A)])
(Array A)
]{
rearranges a Fourier transform @${X} by shifting the zero-frequency component to
the center of the array.

Elements of @${X} are related to coefficients @${c_k} of the trigonometric
interpolant function @${\widetilde{f}} by the formula:

@$${ c_k = \frac{1}{n + 1} \sum_{j=0}^{n}{f(x_j) e^{-ikjh} } }

for @${k = -M, ..., M}, where

@$${ M = \frac{n - 1}{2} }

After @racket[fftshift], the elements of @${X} are rearranged as

@$${ X_{shifted} = [x_{−(M+μ)}, . . . , x_{−1}, x_0, . . . , x_M ] }

}