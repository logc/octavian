#lang racket
(require math/matrix
         math/array)
(provide polyfit
         polyval
         chebyshev-nodes)

;; POLYFIT returns the Lagrange interpolation polynomial of grade `n` that fits
;;   values `ys` at points `xs` with least squared error.
(define (polyfit xs ys grade)
  (define m (add1 grade))
  (define n (length xs))
  (define M
    (build-matrix
      m m
      (λ (row col)
        (for/sum ([i (in-range n)])
          (expt (list-ref xs i) (+ row col))))))
  (define B
    (build-matrix
      m 1
      (λ (row col)
        (for/sum ([i (in-range n)])
          (* (list-ref ys i) (expt (list-ref xs i) row))))))
  (matrix-solve M B))

;; POLYVAL evaluates a polynomial `p` at coordinate `x`.
(define (polyval p x)
  (for/sum ([coeff (in-array p)]
            [exponent (in-range (matrix-num-rows p))])
    (* coeff (expt x exponent))))

;; CHEBYSHEV-NODES are a distribution of n nodes on an interval [a, b] that
;;   yield an optimal Lagrange interpolation polynomial. They are the abscissas
;;   of equispaced nodes on the unit semi-circumference. The Langrage
;;   polynomial resulting from this Chebyshev sampling is less sensitive to
;;   Runge phenomena.
(define (chebyshev-nodes a b n)
  (define x_c (for/list ([i (in-range 0 (add1 n))]) (-1 . * . (cos ((* pi i) . / . n)))))
  (for/list ([x x_c]) ((/ (+ a b) 2) . - . (* (/ (b . - . a) 2) x))))
