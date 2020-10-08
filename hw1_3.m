% 2-norm cond number is cond(A) (or cond(A, norm#))
B = [1 2; 3 4];
C = [1 1; 1 3];
condB = cond(B)
condC = cond(C)
condS = cond(B+C)
condS > condB+condC