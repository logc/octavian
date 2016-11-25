#lang racket

(require rackunit)
(require "nonlinear.rkt")
(test-case
  "Basic bisection tests"
  (let ([𝛆 1e-5]
        [f (lambda (x) x)]
        [g (lambda (x) (+ x 1))])
    (check-= (octavian/bisection f  0 1 𝛆) 0 𝛆)
    (check-= (octavian/bisection f -1 0 𝛆) 0 𝛆)
    (check-= (octavian/bisection f -1 1 𝛆) 0 𝛆)
    (check-= (octavian/bisection g -1 1 𝛆) -1 𝛆)
    (check-= (octavian/bisection g -2 2 𝛆) -1 𝛆)))
(test-case
  "Bisection with max-iter less than correct"
  (let ([𝛆 1e-5]
        [h (lambda (x) (- (sqr x) 1))])
    (check-= (octavian/bisection h -0.25 1.25 𝛆 3) 0.96875 𝛆)
    (check-= (octavian/bisection h -0.25 1.25 𝛆 (octavian/k_min -0.25 1.25 𝛆)) 1 𝛆)))

(test-case
  "Investment fund"
  (let* ([𝛆 (expt 10 -12)] [a 0.01] [b 0.1]
                           [avg-return-rate
                             (lambda (v M r)
                               (M . - . (* v ((1 . + . r) . / . r)
                                           ((expt (1 . + .  r) 5) . - . 1))))]
                           [f (curry avg-return-rate 1000 6000)])
    (check-= (octavian/bisection f a b 𝛆) 0.06140241153618 𝛆)))

(test-case
  "Basic Newton tests"
  (let ([𝛆 1e-5]
        [f (lambda (x) x)]
        [df (lambda (x) 1)]
        [g (lambda (x) (+ x 1))]
        [dg (lambda (x) 1)]
        [h (lambda (x) (- (sqr x) 1))]
        [dh (lambda (x) (* 2 x))])
    (check-= (octavian/newton f df -1 𝛆 1) 0 𝛆)
    (check-= (octavian/newton f df -10 𝛆 1) 0 𝛆)
    (check-= (octavian/newton g dg 0 𝛆 1) -1 𝛆)
    (check-= (octavian/newton g dg 0.22 𝛆 1) -1 𝛆)
    (check-= (octavian/newton h dh 2 𝛆 10) 1 𝛆)
    (check-= (octavian/newton h dh -10 𝛆 10) -1 𝛆)))

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
    (check-= (octavian/newton f df 0.01 𝛆 3) 0.06140241153618 𝛆)))

(test-case
  "Find root 𝜶 = 1 of function `f(x) = e^x * (x - 1)` using Aitken's method"
  (let ([𝜙_0 (lambda (x) (log (* x (exp x))))]
        [𝜙_1 (lambda (x) (/ (+ (exp x) x) (+ (exp x) 1)))]
        [𝛆 1e-10]
        [x_0 2])
    (check-= (octavian/aitken 𝜙_0 x_0 𝛆) 1.0 𝛆)
    (check-= (octavian/aitken 𝜙_1 x_0 𝛆) 1.0 𝛆)))

(test-case
  "Check that Hörner synthetic division does the same as a normal evaluation"
  (check-eqv? (octavian/hoerner '(1 1 1) 2) (+ 1 (* 1 2) (* 1 (sqr 2))))
  (check-eqv? (octavian/hoerner '(0 2 3 1) 3) (+ 0 (* 2 3) (* 3 (sqr 3)) (expt 3 3)))
  (check-eqv? (octavian/hoerner '(1 1 1) 0) 1)
  (check-eqv? (octavian/hoerner '(7 0 1) 1) (+ 7 (sqr 1))))

