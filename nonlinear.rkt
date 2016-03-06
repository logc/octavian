#lang racket

(module+ test
  (require rackunit))

;; BISECTION find a zero of function `f` in interval `[a, b]`, with a tolerance
;;   `𝛆` in the result.  Optionally, pass a maximum number of iterations
;;   `nmax`.  There is a formula to obtain a number of steps `k_min` that are
;;   enough to attain the zero.
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
           (bisection f x b 𝛆 (sub1 nmax))]))]))


;; K_MIN minimum number of steps that are necessary to find a zero in the
;;   interval `[a, b]` with tolerance `𝛆`
(define (k_min a b 𝛆)
  (exact-floor ((logarithm 2 ((- b a) . / . 𝛆)) . - . 1)))


;; LOGARITHM of number `n` in base `base`
(define (logarithm base n)
  (/ (log n) (log base)))


(module+ test

  (test-case
    "Basic bisection tests"
    (let ([𝛆 1e-5]
          [f (lambda (x) x)]
          [g (lambda (x) (+ x 1))])
      (check-= (bisection f  0 1 𝛆) 0 𝛆)
      (check-= (bisection f -1 0 𝛆) 0 𝛆)
      (check-= (bisection f -1 1 𝛆) 0 𝛆)
      (check-= (bisection g -1 1 𝛆) -1 𝛆)
      (check-= (bisection g -2 2 𝛆) -1 𝛆)))


  (test-case
    "Bisection with max-iter less than correct"
    (let ([𝛆 1e-5]
          [h (lambda (x) (- (sqr x) 1))])
      (check-= (bisection h -0.25 1.25 𝛆 3) 0.96875 𝛆)
      (check-= (bisection h -0.25 1.25 𝛆 (k_min -0.25 1.25 𝛆)) 1 𝛆)))


  (test-case
    "Investment fund"
    (let* ([𝛆 (expt 10 -12)] [a 0.01] [b 0.1]
           [avg-return-rate
             (lambda (v M r)
               (M . - . (* v ((1 . + . r) . / . r)
                           ((expt (1 . + .  r) 5) . - . 1))))]
           [f (curry avg-return-rate 1000 6000)])
      (check-= (bisection f a b 𝛆) 0.06140241153618 𝛆))))

;; NEWTON find zeroes of function `f` using Newton's method, which requires the
;;   derivative function `df` and a starting point `x_0`.  The function also
;;   takes a tolerance `𝛆` and a maximum number of iterations `nmax`.
(define (newton f df x_0 𝛆 nmax)
  (let* ([diff (/ (f x_0) (df x_0))]
         [x (- x_0 diff)])
    (cond [(or (= (f x) 0) (= nmax 0) ((magnitude diff) . <= . 𝛆)) x]
          [else (newton f df x 𝛆 (sub1 nmax))])))

(module+ test
  (test-case
    "Basic Newton tests"
    (let ([𝛆 1e-5]
          [f (lambda (x) x)]
          [df (lambda (x) 1)]
          [g (lambda (x) (+ x 1))]
          [dg (lambda (x) 1)]
          [h (lambda (x) (- (sqr x) 1))]
          [dh (lambda (x) (* 2 x))])
      (check-= (newton f df -1 𝛆 1) 0 𝛆)
      (check-= (newton f df -10 𝛆 1) 0 𝛆)
      (check-= (newton g dg 0 𝛆 1) -1 𝛆)
      (check-= (newton g dg 0.22 𝛆 1) -1 𝛆)
      (check-= (newton h dh 2 𝛆 10) 1 𝛆)
      (check-= (newton h dh -10 𝛆 10) -1 𝛆)))

  (test-case
    "Investment fund solved by Newton"
    (let* ([𝛆 1e-12]
          [avg-return-rate (lambda (v M r)
               (M . - . (* v ((1 . + . r) . / . r)
                           ((expt (1 . + .  r) 5) . - . 1))))]
          [davg-return-rate (lambda (v r)
                (+ (* -5 v ((expt (+ r 1) 5) . / . r))
                   (* -1 v (((expt (+ r 1) 5) . - . 1) . / . r))
                   (* v (+ r 1) (((expt (+ r 1) 5) . - . 1) . / .  (sqr r)))))]
          [f (curry avg-return-rate 1000 6000)]
          [df (curry davg-return-rate 1000)])
      (check-= (newton f df 0.01 𝛆 3) 0.06140241153618 𝛆))))

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

(module+ test
  (test-case
    "Find root 𝜶 = 1 of function `f(x) = e^x * (x - 1)` using Aitken's method"
    (let ([𝜙_0 (lambda (x) (log (* x (exp x))))]
          [𝜙_1 (lambda (x) (/ (+ (exp x) x) (+ (exp x) 1)))]
          [𝛆 1e-10]
          [x_0 2])
      (check-= (aitken 𝜙_0 x_0 𝛆) 1.0 𝛆)
      (check-= (aitken 𝜙_1 x_0 𝛆) 1.0 𝛆))))

;; HOERNER evaluates efficiently polynomial `a` at point `z`. The polynomial
;;   is expected to be a list of coefficients (a_0 ... a_n), ordered from lowest
;;   to highest degree.
(define (hoerner a z)
  (cond
    [(empty? a) 0]
    [else (+ (first a) (* z (hoerner (rest a) z)))]))

(module+ test
  (test-case
    "Check that Hörner synthetic division does the same as a normal evaluation"
    (check-eqv? (hoerner '(1 1 1) 2) (+ 1 (* 1 2) (* 1 (sqr 2))))
    (check-eqv? (hoerner '(0 2 3 1) 3) (+ 0 (* 2 3) (* 3 (sqr 3)) (expt 3 3)))
    (check-eqv? (hoerner '(1 1 1) 0) 1)
    (check-eqv? (hoerner '(7 0 1) 1) (+ 7 (sqr 1)))))

;; NEWTONHOERNER computes all roots of the polynomial `a` using Newton-Hoerner's
;;   method, starting at `x_0`. For each root, the method stops after `nmax`
;;   iterations (default: 100) or after the absolute value of the difference
;;   between two consecutive iterations is less than a given tolerance `𝛆`
;;   (default: 1 * 10e-4).
(define (newton-hoerner a x_0 [𝛆 1e-4] [niter 100])
  (cond
    [(= niter 0) (error "Not converged")]
    [(= 1 (length a)) empty]
    [else
      (let ([root (newton (lambda (z) (hoerner a z))
                    (lambda (z) (hoerner (diff a) z))
                    x_0 𝛆 niter)])
        (cons root (newton-hoerner (deflate a root) x_0 𝛆 niter)))]))

(define (diff a [n (length a)])
  (cond
    [(empty? a) empty]
    [(= n (length a)) (diff (rest a) n)]
    [else (cons (* (- n (length a)) (first a))
                (diff (rest a) n))]))

(define (deflate a z)
  (define n (length a))
  (define b (make-list n 0))
  (set! b (list-set b 0 (list-ref a 0)))
  (for ([j (in-range 1 n)])
    (set! b (list-set b j (+ (list-ref a j) (* (list-ref b (sub1 j)) z)))))
  ; return all except the first one, because we express polynomials in reverse
  ; to the reference book
  (rest b))

(module+ test
  ;(test-case
  ;  "Find roots of `p(x) = x^3 - 6x^2 + 11x - 6` using Newton-Hörner method"
  ;  (check-= (first  (newton-hoerner '(-6 11 -6 1) 0 1e-15)) 1 1e-15)
  ;  (check-= (second (newton-hoerner '(-6 11 -6 1) 0 1e-15)) 2 1e-15)
  ;  (check-= (third  (newton-hoerner '(-6 11 -6 1) 0 1e-15)) 3 1e-15)
  ;  )
  (test-case
    "Find roots of `p(x) = x^4 - 7 x^3 + 15 x^2 - 13 x + 4`"
    (check-= (first (newton-hoerner '(4 -13 15 -7 1) 0 1e-15)) 1 1e-5)
    )
  )
