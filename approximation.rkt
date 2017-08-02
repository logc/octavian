#lang racket

(provide polyfit
         polyval
         chebyshev-nodes
         barycentric
         dft
         idft
         cooley-turkey-dft
         rader-dft)

(require math/base
         math/number-theory
         math/matrix
         math/array)

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

;; DFT discrete Fourier transform, not restricted to arrays with a length that
;;   is a power of 2.
(define (dft X)
  (define N (array-size X))
  (define (naive-dft X)
    (for/array ([k (in-range N)])
               (for/sum ([l (in-range N)])
                 (* (array-ref X (vector l))
                    (twiddle-factor N l k)))))
  (cond [(N . < . 100) (naive-dft X)]
        [(power-of-two? N) (array-fft X)]
        [(not  (prime? N)) (cooley-turkey-dft X)]
        [else              (rader-dft X)]))

;; IDFT inverse discrete Fourier transform
(define (idft X)
  (define N (array-size X))
  (define (naive-idft X)
    (for/array ([k (in-range N)])
               (for/sum ([l (in-range N)])
                 (* (array-ref X (vector l))
                    (inverse-twiddle-factor N l k)))))
  (define (scale result)
    (array-map (λ (x) (* (1 . / . N) x)) result))
  (cond [(N . < . 100)     (scale (naive-idft X))]
        [(power-of-two? N) (array-inverse-fft X)]
        [(not  (prime? N)) (scale (cooley-turkey-idft X))]
        [else              (scale (rader-idft X))]))

;; COOLEY-TURKEY-DFT is a fast Fourier transform for arrays of non-prime length
;;   cf. http://cnx.org/contents/ulXtQbN7@15/Implementing-FFTs-in-Practice#uid37
(define (cooley-turkey-dft X)
  (cooley-turkey-algorithm X twiddle-factor))

(define (cooley-turkey-idft X)
  (cooley-turkey-algorithm X inverse-twiddle-factor))

(define (cooley-turkey-algorithm X twiddle-function)
  (define n (array-size X))
  (define-values (n1 n2) (fft-factorize n))
  (for*/array
   ([k2 (in-range n2)]
    [k1 (in-range n1)])
   (for/sum ([l2 (in-range n2)])
     (* (* (for/sum ([l1 (in-range n1)])
             (* (array-ref X (vector (+ (* l1 n2) l2)))
                (twiddle-function n1 l1 k1)))
           (twiddle-function n l2 k1))
        (twiddle-function n2 l2 k2)))))

;; TWIDDLE-FACTOR is a primitive root of unity. In the literature it is often
;;   represented as `w_n^{lk}`
(define (twiddle-factor n l k)
  (ω -1 n l k))

(define (inverse-twiddle-factor n l k)
  (ω +1 n l k))

(define (ω sign n l k)
  (exp (/ (* sign 2 pi 0+i l k) n)))

;; FFT-FACTORIZE returns a small divisor and a big divisor of N. The first one
;;   can be the radix and the second one can be the outer factor for a
;;   Cooley-Turkey FFT algorithm
(define (fft-factorize n)
  (let ([divisor-list (non-trivial-divisors n)])
    (values (first divisor-list) (last divisor-list))))

;; NON-TRIVIAL-DIVISORS returns a list of divisors, excluding the first one, 1,
;;   and the last one, n
(define (non-trivial-divisors n)
  (reverse (cdr (reverse (cdr (divisors n))))))

;; RADER-DFT is a fast Fourier transform for arrays of prime size
;; cf. https://skybluetrades.net/blog/posts/2013/12/31/data-analysis-fft-9.html
;; cf. https://en.wikipedia.org/wiki/Rader%27s_FFT_algorithm
(define (rader-dft X)
  (rader-algorithm X -1))

(define (rader-idft X)
  (rader-algorithm X 1))

(define (rader-algorithm X twiddle-factor-sign)
  (define p (array-size X))
  (define g (primitive-root p))
  ;; index-scheme is a reordering of elements given by g^p in the multiplicative
  ;; group Z_p: https://en.wikipedia.org/wiki/Multiplicative_group_of_integers_modulo_n
  (define index-scheme
    (map (λ (x) (with-modulus p (modexpt g x))) (range (- p 1))))
  (define omega-index-scheme
    (map (λ (x) (with-modulus p (mod/ 1 (modexpt g x)))) (range (- p 1))))
  (define indexes
    (for/array ([idx index-scheme]) (vector idx)))
  (define omega-indexes
    (for/array ([idx omega-index-scheme]) (vector idx)))
  (define inverse-idxs
    (for/list ([k (in-range 1 (add1 (length omega-index-scheme)))]
               [idx omega-index-scheme])
      (cons idx k)))
  (define (find-old-idx new-idx) (cdr (assoc new-idx inverse-idxs)))
  (define (invert-reorder reordered)
    (for/array ([k (in-range 1 (add1 (array-size reordered)))])
               (array-ref reordered (vector (sub1 (find-old-idx k))))))
  (define h_tilde (array-indexes-ref X indexes))
  (define (w x) (ω twiddle-factor-sign p x 1))
  (define omega (list->array (map (lambda (x) (w x)) (range p))))
  (define w_tilde (array-indexes-ref omega omega-indexes))
  (define h_0 (for/sum ([x (in-array X)]) x))
  (define H_tilde (idft (array* (dft h_tilde) (dft w_tilde))))
  (define H (invert-reorder H_tilde))
  (list->array (cons h_0 (array->list H))))
