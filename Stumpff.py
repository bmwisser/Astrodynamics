def Stumpff(z, n_terms):

    def factorial(n):
        if n <= 0:
            return 1
        else:
            return n * factorial(n - 1)

    S = 0
    C = 0

    for i in range(0,n_terms):

        S = S + ((-z)**i/factorial(2*i + 3))
        C = C + ((-z)**i/factorial(2*i + 2))
        #print(S); print(C); print(i); print('\n')
    return S, C

