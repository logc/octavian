#lang racket/base
(provide linspace)

(define (linspace a b n)
  (let ([step (/ (- b a) (exact->inexact n))])
    (build-list (add1 n) (lambda (x) (+ a (* step x))))))
