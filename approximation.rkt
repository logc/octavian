#lang racket
(require math/matrix
         math/array)
(provide polyfit
         polyval
         chebyshev-nodes
         barycentric
         dft)

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
  (define x_c (for/list ([i (in-range 0 (add1 n))])
                (-1 . * . (cos ((* pi i) . / . n)))))
  (for/list ([x x_c]) ((/ (+ a b) 2) . - . (* (/ (b . - . a) 2) x))))

;; BARTYCENTRIC computes the value y_1 at the abscissae x_1 of the polynomial
;;   interpolating data (x, y) by using the barycentric formula.
(define (barycentric xs ys x_1)
  (define n (length xs))
  (define a (apply min xs))
  (define b (apply max xs))
  ;; When n → ∞, barycentric-weights will grow or decay exponentially at the
  ;; rate (4/(b−a))n, where [a,b] is the interval on which we seek the
  ;; interpolating polynomial.
  (define correction (4 . / . (b . - . a)))
  (define (barycentric-weights k)
    (let ([w_k (expt
                (for/product ([x_j xs] [j (in-naturals)] #:unless (= j k))
                  ((list-ref xs k) . - . x_j))
                -1)])
      (* w_k correction)))
  (define barycentric-numerator
    (for/sum ([k (in-range 0 n)] [x_k xs] [y_k ys])
      (((barycentric-weights k) . / . (x_1 . - . x_k)) . * . y_k)))
  (define barycentric-denominator
    (for/sum ([k (in-range 0 n)] [x_k xs])
      ((barycentric-weights k) . / . (x_1 . - . x_k))))
  (barycentric-numerator . / . barycentric-denominator))

;; DFT discrete Fourier transform
(define (dft X)
  (define N (array-size X))
  (define (naive-dft X)
    (for/array ([k (in-range N)])
               (for/sum ([n (in-range N)])
                 (* (array-ref X (vector n))
                    (exp (* (/ (* -2 pi 0+i) N) n k))))))
  (cond [(N . < . 100) (naive-dft X)]
        ;; placeholder for other implementations
        [else (naive-dft X)])
  )
