#lang racket

(module+ test
  (require rackunit))

;; BISECTION find a zero of function `f` in interval `[a, b]`, with a tolerance
;;   `𝛆` in the result.  Optionally, pass a maximum number of iterations
;;   `nmax`.  There is a formula to obtain a number of steps `k_min` that are
;;   enough to attain the zero.
(define (bisection f a b 𝛆 [nmax (k_min a b 𝛆)])
  (cond
    [((f a) . = . 0) a]
    [((f b) . = . 0) b]
    [else
      (let ([x ((+ a b) . / . 2)]
            [I ((- b a) . / . 2)])
        (cond
          [(or (I . < . 𝛆)
               (nmax . = . 0)
               ((f x) . = . 0))
           x]
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
  ((log n) . / . (log base)))


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
  (let* ([diff (* -1 ((f x_0) . / . (df x_0)))]
         [x (x_0 . + . diff)])
    (cond [(or ((f x) . = . 0)
               (nmax . = . 0)
               ((abs diff) . <= . 𝛆))
           x]
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
      (check-= (newton f df 0.01 𝛆 3) 0.06140241153618 𝛆))
    )
  )
