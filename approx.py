#!/usr/bin/env python3

from pathlib import Path
import sys

def sign(i):
    return "+" if i % 2 == 0 else "-"

n = int(sys.argv[1])

directory = f"results/ogr{n}-approx"
Path(directory).mkdir(parents=True, exist_ok=True)

with open(f"results/ogr{n}-approx/stage1.sing", "w") as stage1:
    stage1.write(f"""ring R = (0,v1,v2),z(1..{n}),ws(1..{n});
option(redSB);
""")
    for k in range(n+1,2*n+4):
        stage1.write(f"poly z({k}) = 0;\n")
    for k in range(1,n+1):
        stage1.write(f"poly cs({k}) = 0{sign(k)}(2*z({k})-v1*z({k+1})+2*v1^2*z({k+2})-(8*v1^3+7*v2)*z({k+3}));\n")
    stage1.write("""ideal I =
  z(1)^2-(z(2)+v1*z(3)+v1^2*(2*z(4)-2*z(1)*z(3))+v1^3*(z(5)+4*z(1)*z(4)+z(2)*z(3))+v2*(-z(5)+4*z(1)*z(4)+z(2)*z(3))),
""")
    for k in range(2,n+1):
        constant_term = f"z({2*k})"
        for i in range(1,k):
            constant_term += f"{sign(i)}2*z({i})*z({2*k-i})"
        v1_term = f"{k}*z({2*k+1})-{2*k-2}*z({1})*z({2*k})"
        for i in range(2,k+1):
            v1_term += f"{sign(i)}{2*(k-i)+1}*z({i})*z({2*k+1-i})"
        v1_2_term = f"{(k**2+k-2)//2}*z({2*k+2})-{k**2-k-2}*z({1})*z({2*k+1})+({k**2-2*k-1})*z({2})*z({2*k})"
        for i in range(3,k+1):
            v1_2_term += f"{sign(i)}{(k-i)**2+2*(k-i)+1}*z({i})*z({2*k+2-i})"
        v1_3_term = f"{(k**3+3*k**2+2*k-24)//6}*z({2*k+3})-({(k**3-k-24)//3})*z({1})*z({2*k+2})+({(2*k**3-3*k**2+k-54)//6})*z({2})*z({2*k+1})-({(2*k**3-9*k**2+25*k-72)//6})*z({3})*z({2*k})"
        for i in range(4,k+2):
            v1_3_term += f"{sign(i)}{(2*(k-i)**3+9*(k-i)**2+25*(k-i)+24)//6}*z({i})*z({2*k+3-i})"
        v2_term = f"{k-2}*z({2*k+3})-({2*k-6})*z({1})*z({2*k+2})+({2*k-8})*z({2})*z({2*k+1})-({2*k-10})*z({3})*z({2*k})"
        for i in range(4,k+2):
            v2_term += f"{sign(i)}{2*(k-i)+3}*z({i})*z({2*k+3-i})"
        stage1.write(f"  z({k})^2{sign(k)}(({constant_term})+v1*({v1_term})+v1^2*({v1_2_term})+v1^3*({v1_3_term})+v2*({v2_term})),\n")
    stage1.write(f"""  0;
I = std(I);
""")

with open(f"results/ogr{n}-approx/stage2.sing", "w") as stage2:
    stage2.write(f"< \"results/ogr{n}-approx/stage1.sing\";\n")
    for mask in range(2**(n-1)):
        i_s = list(filter(lambda i: mask & 2**(i-2) != 0, range(2,n+1)))
        for d1 in range(n*(n+1)//2-sum(i_s)+1):
            if d1 + sum(i_s) >= n*(n+1)//2 - 3:
                x = f"z(1)^{d1}"
                for i in i_s:
                    x += f"*cs({i})"
                stage2.write(f"printf(\"{x}=%s\",reduce({x},I));\n")
    stage2.write("quit;\n")
