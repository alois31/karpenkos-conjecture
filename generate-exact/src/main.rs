// SPDX-FileCopyrightText: Alois Wohlschlager <wohlschlager@math.lmu.de>
// SPDX-License-Identifier: EUPL-1.2

use std::{
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
        "ring R = (0,v1,v2,v3,b(16..{d})),cs(1..{n}),ws(1..{n});
poly b(1) = -1;
poly b(2) = -v1;
poly b(3) = -v1^2;
poly b(4) = -2*v1^3-v2;
poly b(5) = -4*v1^4-3*v1*v2;
poly b(6) = -9*v1^5-9*v1^2*v2;
poly b(7) = -21*v1^6-26*v1^3*v2-2*v2^2;
poly b(8) = -51*v1^7-76*v1^4*v2-15*v1*v2^2-v3;
poly b(9) = -127*v1^8-221*v1^5*v2-70*v1^2*v2^2-5*v1*v3;
poly b(10) = -324*v1^9-644*v1^6*v2-281*v1^3*v2^2-7*v2^3-24*v1^2*v3;
poly b(11) = -841*v1^10-1877*v1^7*v2-1044*v1^4*v2^2-77*v1*v2^3-95*v1^3*v3-6*v2*v3;
poly b(12) = -2213*v1^11-5476*v1^8*v2-3708*v1^5*v2^2-505*v1^2*v2^3-339*v1^4*v3-63*v1*v2*v3;
poly b(13) = -5889*v1^12-15997*v1^9*v2-12775*v1^6*v2^2-2618*v1^3*v2^3-1141*v1^5*v3-28*v2^4-399*v1^2*v2*v3;
poly b(14) = -15821*v1^13-46800*v1^10*v2-43061*v1^7*v2^2-11931*v1^4*v2^3-3710*v1^6*v3-414*v1*v2^4-2001*v1^3*v2*v3-40*v2^2*v3;
poly b(15) = -42851*v1^14-137104*v1^11*v2-142794*v1^8*v2^2-50190*v1^5*v2^3-11794*v1^7*v3-3482*v1^2*v2^4-8776*v1^4*v2*v3-556*v1*v2^2*v3-4*v3^2;
poly cs(0) = 1;",
    )
    .unwrap();
    for i in n + 1..=d {
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
option(redSB);
poly d(1) = 2;
poly d(2) = -v1;
poly d(3) = 2*v1^2;
poly d(4) = -8*v1^3-7*v2;
poly d(5) = 26*v1^4+30*v1*v2;
poly d(6) = -84*v1^5-111*v1^2*v2;
poly d(7) = 300*v1^6+502*v1^3*v2+112*v2^2;
poly d(8) = -1140*v1^7-2299*v1^4*v2-960*v1*v2^2-127*v3;
poly d(9) = 4334*v1^8+9958*v1^5*v2+5414*v1^2*v2^2+766*v1*v3;
poly d(10) = -16692*v1^9-43118*v1^6*v2-29579*v1^3*v2^2-2380*v2^3-3579*v1^2*v3;
poly d(11) = 65744*v1^10+189976*v1^7*v2+161034*v1^4*v2^2+31012*v1*v2^3+17770*v1^3*v3+5616*v2*v3;
poly d(12) = -262400*v1^11-837637*v1^8*v2-838452*v1^5*v2^2-240631*v1^2*v2^3-86487*v1^4*v3-55329*v1*v2*v3;
poly d(13) = 1056540*v1^12+3685550*v1^9*v2+4232750*v1^6*v2^2+1600786*v1^3*v2^3+404198*v1^5*v3+58268*v2^4+363210*v1^2*v2*v3;
poly d(14) = -4292816*v1^13-16254540*v1^10*v2-21110372*v1^7*v2^2-10071369*v1^4*v2^3-1864478*v1^6*v3-1022466*v1*v2^4-2193009*v1^3*v2*v3-212440*v2^2*v3;
poly d(15) = 17587492*v1^14+71867828*v1^11*v2+104219628*v1^8*v2^2+60190566*v1^5*v2^3+8581604*v1^7*v3+10170952*v1^2*v2^4+12667346*v1^4*v2*v3+2972696*v1*v2^2*v3+65024*v3^2;\");",
    )
    .unwrap();
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
I = std(I);\");
quit;",
    )
    .unwrap();

    let mut stage3_writer =
        BufWriter::new(File::create(format!("{directory}/stage3.sing")).unwrap());
    writeln!(stage3_writer, "< \"{directory}/stage2.sing\";").unwrap();
    for mask in 0..1 << (n - 1) {
        let is = (2..=n)
            .filter(|i| mask & (1 << (i - 2)) != 0)
            .collect::<Vec<_>>();
        for d1 in 0..=d - is.iter().sum::<u32>() {
            let x = std::iter::once(format!("z(1)^{d1}"))
                .chain(is.iter().map(|i| format!("cs({i})")))
                .collect::<Vec<_>>()
                .join("*");
            writeln!(stage3_writer, "printf(\"{x}=%s\",reduce({x},I));").unwrap();
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
