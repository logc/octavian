#lang racket
(require math/matrix)
(provide polyfit)

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