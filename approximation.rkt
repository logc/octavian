#lang racket
(require math/matrix)
(provide polyfit)

(define (polyfit xs ys grade)
  (define m (add1 grade))
  (define n (length xs)) ; should match (length ys)
  (define M
    (build-matrix
      m m
      (lambda (row col)
        (for/sum ([i (in-range n)])
          (expt (list-ref xs i) (+ row col))))))
  (define B
    (build-matrix
      m 1
      (lambda (row col)
        (for/sum ([i (in-range n)])
          (* (list-ref ys i) (expt (list-ref xs i) row))))))
  (matrix-solve M B))
