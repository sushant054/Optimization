# Write optimization codes based on "Golden ratio search" algorithm.
# Try the function f(x)=x3−3x2+7
#  over interval [1, 3]

import math
# function
def f(x):
    return x**3 - 3*x**2 + 7

def golden_ratio_search(f, a, b, tol=1e-5):
    #golden ration
    gr = (1 + math.sqrt(5)) / 2  # approx. 1.6180..
    # compute c and d
    # b-a=gr*b-gr*c
    # gr*c=gr*b+a-b    <--from this we calculate c and d   
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    while abs(b - a) > tol:
        if f(c) < f(d):
            #if we have   a c d b   then...
            #eleminate rhs from d & replace c to d and d with b
            b = d
        else:
            # if f(c)> f(d) then...
            a = c

        # recalculate c and d
        c = b - (b - a) / gr
        d = a + (b - a) / gr

    return (b + a) / 2
#intervals--> [1, 3]
a = 1
b = 3
minimum = golden_ratio_search(f, a, b)
print(f"minimum value is:{minimum}")


