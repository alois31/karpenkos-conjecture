// SPDX-FileCopyrightText: Alois Wohlschlager <wohlschlager@math.lmu.de>
// SPDX-License-Identifier: EUPL-1.2

use std::{
    borrow::Cow,
    cmp::Reverse,
    collections::HashSet,
    fs::File,
    io::{BufWriter, Write},
};

fn main() {
    let n = std::env::args()
        .nth(1)
        .expect("missing argument")
        .parse::<u32>()
        .expect("invalid number");

    let d = n * (n + 1) / 2;

    let directory = format!("results/ogr{n}");
    std::fs::create_dir_all(&directory).unwrap();
    let mut stage1_writer =
        BufWriter::new(File::create(format!("{directory}/stage1.sing")).unwrap());

    writeln!(
        stage1_writer,
        "ring T = (0,v1,v2,v3),t,ls;
poly l = t+1/2*v1*t^2+(1/4*v1^3+1/2*v2)*t^4+(1/8*v1^7+1/4*v1^4*v2+1/4*v1*v2^2+1/2*v3)*t^8;
poly e = t-1/2*v1*t^2+1/2*v1^2*t^3+(-7/8*v1^3-1/2*v2)*t^4+(13/8*v1^4+3/2*v1*v2)*t^5+(-49/16*v1^5-7/2*v1^2*v2)*t^6+(97/16*v1^6+17/2*v1^3*v2+v2^2)*t^7+(-1615/128*v1^7-683/32*v1^4*v2-47/8*v1*v2^2-1/2*v3)*t^8+(3457/128*v1^8+1701/32*v1^5*v2+175/8*v1^2*v2^2+5/2*v1*v3)*t^9+(-15059/256*v1^9-2101/16*v1^6*v2-1133/16*v1^3*v2^2-11/4*v2^3-33/4*v1^2*v3)*t^10+(33311/256*v1^10+2599/8*v1^7*v2+1747/8*v1^4*v2^2+97/4*v1*v2^3+97/4*v1^3*v3+3*v2*v3)*t^11+(-298779/1024*v1^11-206401/256*v1^8*v2-20839/32*v1^5*v2^2-1001/8*v1^2*v2^3-273/4*v1^4*v3-91/4*v1*v2*v3)*t^12+(677425/1024*v1^12+513051/256*v1^9*v2+60417/32*v1^6*v2^2+2065/4*v1^3*v2^3+1491/8*v1^5*v3+35/4*v2^4+105*v1^2*v2*v3)*t^13+(-3100153/2048*v1^13-1276797/256*v1^10*v2-343305/64*v1^7*v2^2-30585/16*v1^4*v2^3-7965/16*v1^6*v3-825/8*v1*v2^4-795/2*v1^3*v2*v3-15*v2^2*v3)*t^14+(7150065/2048*v1^14+3181831/256*v1^11*v2+1923287/128*v1^8*v2^2+106147/16*v1^5*v2^3+20977/16*v1^7*v3+5461/8*v1^2*v2^4+5465/4*v1^4*v2*v3+155*v1*v2^2*v3+2*v3^2)*t^15;
proc coeffs_t(poly p) {{
    matrix m = coeffs(p,t);
    matrix n[15][1] = m[2..nrows(m),1];
    return(n);
}}
matrix b = coeffs_t(subst(e,t,-l));
matrix d = coeffs_t(subst(e,t,2*l));
matrix l_ = coeffs_t(l);
matrix e_ = coeffs_t(e);
ring R = (0,v1,v2,v3,b(16..{d})),cs(1..{n}),ws(1..{n});
matrix b = fetch(T,b);
matrix d = fetch(T,d);
matrix l_ = fetch(T,l_);
matrix e_ = fetch(T,e_);",
    ).unwrap();
    for i in 1..=15 {
        writeln!(stage1_writer, "poly b({i}) = b[{i},1];").unwrap();
    }
    writeln!(stage1_writer, "poly cs(0) = 1;").unwrap();
    for i in n + 1..=u32::max(d, 15) {
        writeln!(stage1_writer, "poly cs({i}) = 0;").unwrap();
    }

    // Step 1: compute expressions for the Chern classes of the dual
    let mut known_monomial_symmetric_functions = HashSet::new();
    writeln!(stage1_writer, "poly c(0) = 1;").unwrap();
    for i in 1..=n {
        let ci = (i..=d)
            .flat_map(|w| partitions(w, i))
            .map(|j| {
                j.iter()
                    .map(|part| format!("b({part})"))
                    .chain(std::iter::once(define_monomial_symmetric_function(
                        n,
                        &j,
                        &mut known_monomial_symmetric_functions,
                        &mut stage1_writer,
                    )))
                    .collect::<Vec<_>>()
                    .join("*")
            })
            .collect::<Vec<_>>()
            .join("+");
        writeln!(stage1_writer, "poly c({i}) = {ci};").unwrap();
    }

    // Step 2: compute the Chern subalgebra
    writeln!(stage1_writer, "ideal I =").unwrap();
    for k in (1..=n).rev() {
        writeln!(
            stage1_writer,
            "  {},",
            ((2 * k).saturating_sub(n)..=u32::min(n, 2 * k))
                .map(|i| format!("cs({i})*c({})", 2 * k - i))
                .collect::<Vec<_>>()
                .join("+"),
        )
        .unwrap();
    }
    writeln!(
        stage1_writer,
        "  cs({n})^2;
I = std(I);",
    )
    .unwrap();

    // Step 3: write the next stage, which will compute the cohomology ring
    writeln!(
        stage1_writer,
        "print(\"ring R = (0,v1,v2,v3,b(16..{d}),d(16..{n})),z(1..{n}),ws(1..{n});
option(redSB);\");",
    )
    .unwrap();
    for i in 1..=15 {
        writeln!(stage1_writer, "printf(\"poly d({i}) = %s;\",d[{i},1]);").unwrap();
    }
    for i in 1..=n {
        writeln!(
            stage1_writer,
            "print(\"poly cs({i}) = (-1)^{i}*({});\");",
            (0..=(n - i))
                .map(|k| format!("d({})*z({})", k + 1, k + i))
                .collect::<Vec<_>>()
                .join("+"),
        )
        .unwrap();
    }
    writeln!(stage1_writer, "print(\"ideal I =\");").unwrap();
    for i in 1..=n {
        writeln!(
            stage1_writer,
            "printf(\"  cs({i})^2-(%s),\",reduce(cs({i})^2,I));",
        )
        .unwrap();
    }
    writeln!(
        stage1_writer,
        "print(\"  0;
I = std(I);\");",
    )
    .unwrap();
    for i in 1..=15 {
        writeln!(
            stage1_writer,
            "poly p({i}) = {};",
            (1..i)
                .map(|j| format!("(-1)^{}*cs({j})*p({})", j - 1, i - j))
                .chain(std::iter::once(format!("(-1)^{}*{i}*cs({i})", i - 1)))
                .collect::<Vec<_>>()
                .join("+")
        )
        .unwrap();
    }
    writeln!(
        stage1_writer,
        "poly lu = reduce(({})/2,I);
poly u = 0;
for (int i=15; i>0; i--) {{
    u = reduce(lu*u+e_[i,1],I);
}}
u = reduce(lu*u,I);
printf(\"poly u = reduce(%s,I);\",u);",
        (1..=15)
            .map(|i| format!("l_[{i},1]*p({i})"))
            .collect::<Vec<_>>()
            .join("+"),
    )
    .unwrap();
    writeln!(stage1_writer, "quit;").unwrap();

    let mut stage3_writer =
        BufWriter::new(File::create(format!("{directory}/stage3.sing")).unwrap());
    writeln!(
        stage3_writer,
        "< \"{directory}/stage2.sing\";
poly x;",
    )
    .unwrap();
    for mask in 0..1 << (n - 1) {
        let is = (2..=n)
            .filter(|i| mask & (1 << (i - 2)) != 0)
            .collect::<Vec<_>>();
        for d1 in 0..=d - is.iter().sum::<u32>() {
            writeln!(stage3_writer, "x = 1;").unwrap();
            for factor in std::iter::repeat_n(Cow::Borrowed("u"), d1.try_into().unwrap())
                .chain(is.iter().map(|i| Cow::Owned(format!("cs({i})"))))
                .rev()
            {
                writeln!(stage3_writer, "x = reduce({factor}*x,I);").unwrap();
            }
            writeln!(
                stage3_writer,
                "printf(\"{}=%s\",x);",
                std::iter::once(format!("u^{d1}"))
                    .chain(is.iter().map(|i| format!("cs({i})")))
                    .collect::<Vec<_>>()
                    .join("*"),
            )
            .unwrap();
        }
    }
    writeln!(stage3_writer, "quit;").unwrap();
}

fn partitions(weight: u32, length: u32) -> Vec<Vec<u32>> {
    if length == 1 {
        vec![vec![weight]]
    } else {
        (1..=weight - length + 1)
            .flat_map(|i| {
                partitions(weight - i, length - 1)
                    .into_iter()
                    .filter_map(move |j| (j[0] <= i).then(|| std::iter::once(i).chain(j).collect()))
            })
            .collect()
    }
}

fn define_monomial_symmetric_function(
    n: u32,
    j: &[u32],
    known: &mut HashSet<Vec<u32>>,
    writer: &mut dyn Write,
) -> String {
    let name = format!(
        "m({})",
        j.iter()
            .map(|part| format!("{part}"))
            .collect::<Vec<_>>()
            .join(",")
    );
    if !known.contains(j) {
        if j.len() > n as usize {
            writeln!(writer, "poly {name} = 0;").unwrap();
        } else if j[0] == 1 {
            writeln!(writer, "poly {name} = cs({});", j.len()).unwrap();
        } else {
            let jr = j
                .iter()
                .filter_map(|part| (*part > 1).then_some(part - 1))
                .collect::<Vec<_>>();
            let mjr = define_monomial_symmetric_function(n, &jr, known, writer);
            let m_neighbours = neighbours(&jr, j.len())
                .into_iter()
                .filter(|jn| jn != j)
                .map(|jn| {
                    format!(
                        "{}*{}",
                        neighbour_multiplicity(&jn, &jr),
                        define_monomial_symmetric_function(n, &jn, known, writer),
                    )
                })
                .collect::<Vec<_>>()
                .join("+");
            writeln!(
                writer,
                "poly {name} = reduce(cs({})*{mjr}-({m_neighbours}),cs({n})^2);",
                j.len(),
            )
            .unwrap();
        }
        known.insert(j.to_owned());
    }
    name
}

fn neighbours(jr: &[u32], k: usize) -> HashSet<Vec<u32>> {
    if jr.is_empty() {
        HashSet::from_iter([vec![1; k]])
    } else {
        neighbours(&jr[1..], k - 1)
            .into_iter()
            .map(|jn| std::iter::once(jr[0] + 1).chain(jn).collect())
            .chain(neighbours(&jr[1..], k).into_iter().map(|jn| {
                let mut result = std::iter::once(jr[0]).chain(jn).collect::<Vec<_>>();
                result.sort_by_key(|part| Reverse(*part));
                result
            }))
            .collect()
    }
}

fn neighbour_multiplicity(jn: &[u32], jr: &[u32]) -> u32 {
    do_neighbour_multiplicity(
        jn,
        &jr.iter()
            .copied()
            .chain(std::iter::repeat_n(0, jn.len() - jr.len()))
            .collect::<Vec<_>>(),
    )
}

fn do_neighbour_multiplicity(jn: &[u32], jr: &[u32]) -> u32 {
    if jn.is_empty() {
        1
    } else {
        assert!(jn[0] == jr[0] || jn[0] == jr[0] + 1);
        let lead_n = jn.iter().take_while(|part| **part == jn[0]).count();
        let lead_r = jr.iter().take_while(|part| **part == jr[0]).count();
        if jn[0] > jr[0] {
            assert!(lead_n <= lead_r);
            neighbour_multiplicity(&jn[lead_n..], &jr[lead_n..])
        } else {
            assert!(lead_n >= lead_r);
            assert!((lead_r..lead_n).all(|i| jn[i] == jr[i] + 1));
            binom(lead_n as u32, lead_r as u32)
                * neighbour_multiplicity(&jn[lead_n..], &jr[lead_n..])
        }
    }
}

fn binom(n: u32, k: u32) -> u32 {
    if k == 0 {
        1
    } else {
        n * binom(n - 1, k - 1) / k
    }
}
