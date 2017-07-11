#lang racket
(provide (all-defined-out))

;; BISECTION find a zero of function `f` in interval `[a, b]`, with a tolerance
;;   `𝛆` in the result.  Optionally, pass a maximum number of iterations
;;   `nmax`.  There is a formula to obtain a number of steps `k_min` that are
;;   enough to reach the zero.
(define (bisection f a b 𝛆 [nmax (k_min a b 𝛆)])
  (cond
    [(= (f a) 0) a]
    [(= (f b) 0) b]
    [else
      (let ([x (/ (+ a b) 2)]
            [I (/ (- b a) 2)])
        (cond
          [(or (I . < . 𝛆) (= nmax 0) (= (f x) 0)) x]
          [((* (f a) (f x)) . < . 0)
           (bisection f a x 𝛆 (sub1 nmax))]
          [((* (f x) (f b)) . < . 0)
           (bisection f x b 𝛆 (sub1 nmax))]
          [else error "could not find a zero"]))]))

;; K_MIN minimum number of steps that are necessary to find a zero in the
;;   interval `[a, b]` with tolerance `𝛆`
(define (k_min a b 𝛆)
  (exact-floor ((logarithm 2 ((- b a) . / . 𝛆)) . - . 1)))

;; LOGARITHM of number `n` in base `base`
(define (logarithm base n)
  (/ (log n) (log base)))


;; NEWTON find zeroes of function `f` using Newton's method, which requires the
;;   derivative function `df` and a starting point `x_0`.  The function also
;;   takes a tolerance `𝛆` and a maximum number of iterations `nmax`.
(define (newton f df x_0 𝛆 nmax)
  (let* ([diff (/ (f x_0) (df x_0))]
         [x (- x_0 diff)])
    (cond [(or (= (f x) 0) (= nmax 0) ((magnitude diff) . <= . 𝛆)) x]
          [else (newton f df x 𝛆 (sub1 nmax))])))

;; NEWTONSYS not implemented
(define (newtonsys) empty)

;; AITKEN finds approximations of fixed point `𝜶` of function `𝜙` starting from
;;   initial point `x_0` using Aitken's extrapolation method. The method stops
;;   after a given `nmax` iterations (default: 100), or after the absolute
;;   value of the difference between two consecutive iterations is less than a
;;   given tolerance `𝛆` (default: 1 * 10e-4).
(define (aitken 𝜙 x [𝛆 1e-4] [nmax 100] [niter 0] [diff (add1 𝛆)])
  (cond
    [(or (>= niter nmax) (<= diff 𝛆)) x]
    [else (let* ([gx (𝜙 x)]
                 [ggx (𝜙 gx)]
                 [x_new (/ (- (* x ggx) (sqr gx)) (+ ggx (* -2 gx) x))])
            (aitken 𝜙 x_new 𝛆 nmax (add1 niter) (magnitude (- x x_new))))]))

;; HOERNER evaluates efficiently polynomial `a` at point `z`. The polynomial
;;   is expected to be a list of coefficients (a_0 ... a_n), ordered from lowest
;;   to highest degree.
(define (hoerner a z)
  (cond
    [(empty? a) 0]
    [else (+ (first a) (* z (hoerner (rest a) z)))]))
