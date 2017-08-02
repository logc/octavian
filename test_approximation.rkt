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

(define (array-contiguous-slice arr start end)
  (array-indexes-ref arr (for/array ([k (in-range start end)]) (vector k))))

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
  (check-array-= (dft (array #[1 1])) (array #[2 0]) 1e-9)
  (let* ([n 9]
         [xs (for/list ([k (in-range (+ n 1))]) (* 2 pi (1 . / . (+ n 1)) k))]
         [ys (map (lambda (x) (* x (x . - . (* 2 pi)) (exp (* -1 x)))) xs)]
         [Y (dft (list->array ys))]
         [Y-ct (cooley-turkey-dft (list->array ys))])
    (check-array-= Y-ct Y 1e-9))
  (let* ([n 10]
         [xs (for/list ([k (in-range (+ n 1))]) (* 2 pi (1 . / . (+ n 1)) k))]
         [ys (map (lambda (x) (* x (x . - . (* 2 pi)) (exp (* -1 x)))) xs)]
         [actual (idft (list->array ys))]
         ;; Numbers taken from Octave's ifft function
         [expected (array #[-0.65749
                            -0.05236-0.42048i
                            0.12061-0.16318i
                            0.10256-0.06206i
                            0.08339-0.02501i
                            0.07455-0.00690i
                            0.07455+0.00690i
                            0.08339+0.02501i
                            0.10256+0.06206i
                            0.12061+0.16318i
                            -0.05236+0.42048i])])
    (check-array-= actual expected 1e-5))
  (let* ([n 10]
         [xs (for/list ([k (in-range (+ n 1))]) (* 2 pi (1 . / . (+ n 1)) k))]
         [ys (map (lambda (x) (* x (x . - . (* 2 pi)) (exp (* -1 x)))) xs)]
         [actual (rader-dft (list->array ys))]
         ;; Numbers taken from Octave's fft function
         [expected (array #[-7.23237
                            -0.57596+4.62532i
                            1.32667+1.79495i
                            1.12813+0.68268i
                            0.91729+0.27506i
                            0.82007+0.07585i
                            0.82007-0.07585i
                            0.91729-0.27506i
                            1.12813-0.68268i
                            1.32667-1.79495i
                            -0.57596-4.62532i])])
    (check-array-= actual expected 1e-5))
  ;; Cooley-Turkey algorithm
  (let* ([n 104]
         [xs (for/list ([k (in-range (+ n 1))]) (* 2 pi (1 . / . (+ n 1)) k))]
         [ys (map (lambda (x) (* x (x . - . (* 2 pi)) (exp (* -1 x)))) xs)]
         [actual (dft (list->array ys))]
         [expected (array #[-71.8046
                            -8.3086+44.2580i
                            9.7192+17.3651i
                            7.5798+6.9122i
                            5.1723+3.2650i
                            ])])
    (check-array-= (array-contiguous-slice actual 0 5) expected 1e-4))
  ;; Radix-2 DIT implemented in array-fft
  (let* ([n 127]
         [xs (for/list ([k (in-range (+ n 1))]) (* 2 pi (1 . / . (+ n 1)) k))]
         [ys (map (lambda (x) (* x (x . - . (* 2 pi)) (exp (* -1 x)))) xs)]
         [actual (dft (list->array ys))]
         [expected (array #[-87.5457
                            -10.1411+53.9526i
                            11.8357+21.1689i
                            9.2275+8.4263i
                            6.2928+3.9803i])])
    (check-array-= (array-contiguous-slice actual 0 5) expected 1e-4))
  ;; Rader algorithm
  (let* ([n 100]
         [xs (for/list ([k (in-range (+ n 1))]) (* 2 pi (1 . / . (+ n 1)) k))]
         [ys (map (lambda (x) (* x (x . - . (* 2 pi)) (exp (* -1 x)))) xs)]
         [actual (dft (list->array ys))]
         [expected (array #[-69.0667
                            -7.9897+42.5720i
                            9.3514+16.7036i
                            7.2935+6.6489i
                            4.9778+3.1406i])])
    (check-array-= (array-contiguous-slice actual 0 5) expected 1e-4))
  (let* ([n 100]
         [xs (for/list ([k (in-range (+ n 1))]) (* 2 pi (1 . / . (+ n 1)) k))]
         [ys (map (lambda (x) (* x (x . - . (* 2 pi)) (exp (* -1 x)))) xs)]
         [actual (idft (list->array ys))]
         [expected (array #[-0.68382
                            -0.07910-0.42150i
                            0.09259-0.16538i
                            0.07222-0.06583i
                            0.04929-0.03110i
                            ])])
    (check-array-= (array-contiguous-slice actual 0 5) expected 1e-4))
  )


