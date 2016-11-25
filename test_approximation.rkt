#lang racket

(require rackunit)
(require math/matrix)
(require "approximation.rkt")

(define-simple-check (check-all-= actual-numbers expected-numbers epsilon)
                     (for ([actual actual-numbers]
                           [expected expected-numbers])
                       (check-= actual expected epsilon)))

(define-simple-check (check-matrix-= M N epsilon)
                     (check-all-= (matrix->list M) (matrix->list N) epsilon))

(test-case
  "Climatology interpolation example"
  (let ([xs '(-55 -25 5 35 65)]
        [ys '(-3.25 -3.2 -3.02 -3.32 -3.1)])
    (check-matrix-= (polyfit xs ys 4)
                    (col-matrix [-3.0132 3.7757e-4 -3.4684e-4 -4.5267e-7 8.2819e-8])
                    1e-3)))
