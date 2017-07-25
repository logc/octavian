#lang racket

(require rackunit)
(require math/matrix
         math/array)
(require "approximation.rkt"
         "fundamentals.rkt")

(define-simple-check (check-all-= actual-numbers expected-numbers epsilon)
  (for ([actual actual-numbers]
        [expected expected-numbers])
    (check-= actual expected epsilon)))

(define-simple-check (check-matrix-= M N epsilon)
  (check-all-= (matrix->list M) (matrix->list N) epsilon))

(define-simple-check (check-array-= A B epsilon)
  (check-all-= (array->list A) (array->list B) epsilon))

(test-case
    "Polyval"
  (define p_0 (col-matrix [1.0]))
  (define p_1 (col-matrix [0.0 2.0]))
  (define p_11 (col-matrix [1.0 3.0]))
  (define p_2 (col-matrix [0.0 0.0 1.0]))
  (define p_22 (col-matrix [1.0 0.0 1.0]))
  (define p_23 (col-matrix [2.0 1.0 3.0]))
  (check-eqv? (polyval p_0 1.0) 1.0)
  (check-eqv? (polyval p_0 (random 1 100)) 1.0)
  (check-eqv? (polyval p_1 2.0) 4.0)
  (check-eqv? (polyval p_11 2.0) (+ 1.0 (* 3.0 2.0)))
  (check-eqv? (polyval p_2 2.0) 4.0)
  (check-eqv? (polyval p_2 3.0) 9.0)
  (check-eqv? (polyval p_22 2.0) (+ 1.0 (sqr 2.0)))
  (check-eqv? (polyval p_23 3.0) (+ 2.0 3.0 (* 3.0 (sqr 3.0)))))

(test-case
  "Climatology interpolation example"
  (let ([xs '(-55 -25 5 35 65)]
        [ys '(-3.25 -3.2 -3.02 -3.32 -3.1)])
    (check-matrix-= (polyfit xs ys 4)
                    (col-matrix [-3.0132 3.7757e-4 -3.4684e-4 -4.5267e-7 8.2819e-8])
                    1e-3)))

(test-case
    "Chebyshev nodes"
  (let ([a -5]
        [b 5]
        [f (λ (x) (1 . / . (1 . + . (sqr x))))]
        )
    (define (err n)
      (define xs (chebyshev-nodes a b n))
      (define ys (map f xs))
      (define c (polyfit xs ys n))
      (define x1 (linspace -5 5 1000))
      (define p (map (λ (x) (polyval c x)) x1))
      (define f1 (map f x1))
      (apply max (map (compose abs -) p f1)))
    (check-= (err  5) 0.6386 1e-4)
    (check-= (err 10) 0.1322 1e-4)
    (check-= (err 20) 0.0177 1e-4)
    ;; this example outputs the following warning in octave 4.2.1:
    ;; warning: matrix singular to machine precision, rcond = 1.42019e-31
    #;(check-= (err 40) 3.3996e-04 1e-4)))

(test-case
    "Barycentric interpolation"
  (define a -5)
  (define b 5)
  (define f (λ (x) (1 . / . (1 . + . (sqr x)))))
  (let* ([xs (chebyshev-nodes a b 4)]
         [ys (map f xs)])
    (check-= (barycentric xs ys 4.9) 0.0088179 1e-7)
    (check-= (barycentric xs ys -1.0) 0.89316 1e-5)
    (check-= (barycentric xs ys 0.0) 1.0 1e-5)
    (check-= (barycentric xs ys 227) 7.5591e+06 1e2))
  (let* ([xs (chebyshev-nodes a b 12)]
         [ys (map f xs)])
    (check-= (barycentric xs ys 4.9) 0.035205 1e-6)
    (check-= (barycentric xs ys -1.0) 0.55924 1e-5)
    (check-= (barycentric xs ys 0.0) 1.0 1e-5)
    ;; Chosen nearer to b because further away the implementations diverge (?)
    (check-= (barycentric xs ys 22) 3.2619e9 1e5))
  (let* ([xs (chebyshev-nodes a b 128)]
         [ys (map f xs)])
    (check-= (barycentric xs ys 4.9) 0.039984 1e-6)
    (check-= (barycentric xs ys -1.0) 1/2 1e-6)
    (check-= (barycentric xs ys 0.0) 1.0 1e-9)
    ;; Notice how this has to be even nearer to b
    (check-= (barycentric xs ys 5.1) 0.079201 1e-6)))

(test-case
    "Trigonometric interpolation and FFT"
  (check-array-= (dft (array #[1])) (array #[1]) 1e-9)
  (check-array-= (dft (array #[1 1])) (array #[2 0]) 1e-9))
