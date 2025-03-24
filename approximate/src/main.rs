// SPDX-FileCopyrightText: Alois Wohlschlager <wohlschlager@math.lmu.de>
// SPDX-License-Identifier: EUPL-1.2

use rustc_hash::FxHasher;
use std::{
    collections::HashMap,
    fmt::{Display, Formatter},
    hash::BuildHasherDefault,
    ops::BitOr,
};

type Exponent = u8;
type Coefficient = i16;
type ZMask = u32;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct BasisElement {
    v1: Exponent,
    v2: Exponent,
    zs: ZMask,
}

impl Display for BasisElement {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self.v1 {
            0 => {}
            1 => write!(f, "v_1")?,
            e => write!(f, "v_1^{e}")?,
        }
        match self.v2 {
            0 => {}
            1 => write!(f, "v_2")?,
            e => write!(f, "v_2^{e}")?,
        }
        write!(
            f,
            "z_{{{}}}",
            (1..=ZMask::BITS)
                .filter(|i| self.zs & (1 << (i - 1)) != 0)
                .map(|i| i.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

type Element = HashMap<BasisElement, Coefficient, BuildHasherDefault<FxHasher>>;

fn term(
    c: Coefficient,
    v1: Exponent,
    v2: Exponent,
    zs: &[u8],
    n: u8,
) -> Option<(BasisElement, Coefficient)> {
    zs.iter().all(|i| *i <= n).then(|| {
        (
            BasisElement {
                v1,
                v2,
                zs: zs.iter().map(|i| 1 << (i - 1)).fold(0, BitOr::bitor),
            },
            c,
        )
    })
}

fn multiply_monomial_z(
    m: BasisElement,
    k: u8,
    t: u8,
    squares: &HashMap<u8, Element>,
    cache: &mut HashMap<(BasisElement, u8), Element>,
) -> Element {
    if cache.contains_key(&(m, k)) {
        return cache[&(m, k)].clone();
    }
    let result = if m.zs & (1 << (k - 1)) == 0 {
        Element::from_iter(std::iter::once((
            BasisElement {
                v1: m.v1,
                v2: m.v2,
                zs: m.zs | (1 << (k - 1)),
            },
            1,
        )))
    } else {
        multiply_element_monomial(
            squares[&k].clone(),
            BasisElement {
                v1: m.v1,
                v2: m.v2,
                zs: m.zs & !(1 << (k - 1)),
            },
            t,
            squares,
            cache,
        )
    };
    cache.insert((m, k), result.clone());
    result
}

fn multiply_monomial_monomial(
    m1: BasisElement,
    m2: BasisElement,
    t: u8,
    squares: &HashMap<u8, Element>,
    cache: &mut HashMap<(BasisElement, u8), Element>,
) -> Element {
    if m2.zs == 0 {
        let v1 = m1.v1 + m2.v1;
        let v2 = m1.v2 + m2.v2;
        if v1 + 3 * v2 > 3 {
            Element::default()
        } else {
            Element::from_iter(std::iter::once((BasisElement { v1, v2, zs: m1.zs }, 1)))
        }
    } else {
        let k = u8::try_from(m2.zs.ilog2() + 1).unwrap();
        multiply_element_monomial(
            multiply_monomial_z(m1, k, t, squares, cache),
            BasisElement {
                v1: m2.v1,
                v2: m2.v2,
                zs: m2.zs & !(1 << (k - 1)),
            },
            t,
            squares,
            cache,
        )
    }
}

fn multiply_element_monomial(
    e: Element,
    m2: BasisElement,
    t: u8,
    squares: &HashMap<u8, Element>,
    cache: &mut HashMap<(BasisElement, u8), Element>,
) -> Element {
    let mut result = Element::default();
    for (m1, c1) in e {
        for (m, c) in multiply_monomial_monomial(m1, m2, t, squares, cache) {
            let entry = result.entry(m).or_insert(0);
            *entry = (*entry + c1 * c) & ((1 << t) - 1);
        }
    }
    result.retain(|_, c| *c != 0);
    result
}

fn multiply_element_element(
    e1: Element,
    e2: Element,
    t: u8,
    squares: &HashMap<u8, Element>,
    cache: &mut HashMap<(BasisElement, u8), Element>,
) -> Element {
    let mut result = Element::default();
    for (m2, c2) in e2 {
        for (m, c) in multiply_element_monomial(e1.clone(), m2, t, squares, cache) {
            let entry = result.entry(m).or_insert(0);
            *entry = (*entry + c * c2) & ((1 << t) - 1);
        }
    }
    result.retain(|_, c| *c != 0);
    result
}

fn rational_element(
    k: Exponent,
    cs: &[u8],
    n: u8,
    t: u8,
    squares: &HashMap<u8, Element>,
    cache: &mut HashMap<(BasisElement, u8), Element>,
) -> Element {
    std::iter::repeat_n(Element::from_iter(term(1, 0, 0, &[1], n)), k.into())
        .chain(cs.iter().map(|k| {
            [
                term(2, 0, 0, &[*k], n),
                term(-1, 1, 0, &[k + 1], n),
                term(2, 2, 0, &[k + 2], n),
                term(-8, 3, 0, &[k + 3], n),
                term(-7, 0, 1, &[k + 3], n),
            ]
            .into_iter()
            .flatten()
            .map(|(m, c)| {
                (
                    m,
                    (Coefficient::pow(-1, (k + 1).into()) * c) & ((1 << t) - 1),
                )
            })
            .filter(|(m, c)| *c != 0 && m.zs == m.zs & ((1 << n) - 1))
            .collect()
        }))
        .rfold(Element::from_iter(term(1, 0, 0, &[], n)), |accum, e| {
            multiply_element_element(accum, e, t, squares, cache)
        })
}

fn torsion_exponent(n: u8) -> u8 {
    if n == 0 {
        return 0;
    }
    let n_ = u16::from(n);
    let temp = n_ - u16::try_from((n_ * (n_ + 1) / 2 + 1).ilog2()).unwrap();
    let e = u8::try_from(n.ilog2()).unwrap();
    let b = n - (1 << e);
    let exceptional = 2 * b + 3 <= e + torsion_exponent(b);
    u8::try_from(temp).unwrap() + if exceptional { 1 } else { 0 }
}

fn main() {
    let n = std::env::args()
        .nth(1)
        .expect("missing argument")
        .parse()
        .expect("invalid number");

    let t = torsion_exponent(n);
    let squares = std::iter::once((
        1,
        [
            term(1, 0, 0, &[2], n),
            term(1, 1, 0, &[3], n),
            term(2, 2, 0, &[4], n),
            term(-2, 2, 0, &[1, 3], n),
            term(1, 3, 0, &[5], n),
            term(4, 3, 0, &[1, 4], n),
            term(1, 3, 0, &[2, 3], n),
            term(-1, 0, 1, &[5], n),
            term(4, 0, 1, &[1, 4], n),
            term(1, 0, 1, &[2, 3], n),
        ]
        .into_iter()
        .flatten()
        .map(|(m, c)| (m, c & ((1 << t) - 1)))
        .filter(|(m, c)| *c != 0 && m.zs == m.zs & ((1 << n) - 1))
        .collect(),
    ))
    .chain((2..=n).map(|k| {
        let k_ = Coefficient::from(k);
        (
            k,
            term(1, 0, 0, &[2 * k], n)
                .into_iter()
                .chain((1..=k - 1).flat_map(|i| {
                    term(Coefficient::pow(-1, i.into()) * 2, 0, 0, &[i, 2 * k - i], n)
                }))
                .chain(term(k_, 1, 0, &[2 * k + 1], n))
                .chain(term(-(2 * k_ - 2), 1, 0, &[1, 2 * k], n))
                .chain((2..=k).flat_map(|i| {
                    let i_ = Coefficient::from(i);
                    term(
                        Coefficient::pow(-1, i.into()) * (2 * (k_ - i_) + 1),
                        1,
                        0,
                        &[i, 2 * k + 1 - i],
                        n,
                    )
                }))
                .chain(term((k_ * k_ + k_ - 2) / 2, 2, 0, &[2 * k + 2], n))
                .chain(term(-(k_ * k_ - k_ - 2), 2, 0, &[1, 2 * k + 1], n))
                .chain(term(k_ * k_ - 2 * k_ - 1, 2, 0, &[2, 2 * k], n))
                .chain((3..=k).flat_map(|i| {
                    let i_ = Coefficient::from(i);
                    term(
                        Coefficient::pow(-1, i.into())
                            * ((k_ - i_) * (k_ - i_) + 2 * (k_ - i_) + 1),
                        2,
                        0,
                        &[i, 2 * k + 2 - i],
                        n,
                    )
                }))
                .chain(term(
                    (k_ * k_ * k_ + 3 * k_ * k_ + 2 * k_ - 24) / 6,
                    3,
                    0,
                    &[2 * k + 3],
                    n,
                ))
                .chain(term(
                    -(k_ * k_ * k_ - k_ - 24) / 3,
                    3,
                    0,
                    &[1, 2 * k + 2],
                    n,
                ))
                .chain(term(
                    (2 * k_ * k_ * k_ - 3 * k_ * k_ + k_ - 54) / 6,
                    3,
                    0,
                    &[2, 2 * k + 1],
                    n,
                ))
                .chain(term(
                    -(2 * k_ * k_ * k_ - 9 * k_ * k_ + 25 * k_ - 72) / 6,
                    3,
                    0,
                    &[3, 2 * k],
                    n,
                ))
                .chain((4..=k + 1).flat_map(|i| {
                    let i_ = Coefficient::from(i);
                    term(
                        Coefficient::pow(-1, i.into())
                            * (2 * (k_ - i_) * (k_ - i_) * (k_ - i_)
                                + 9 * (k_ - i_) * (k_ - i_)
                                + 25 * (k_ - i_)
                                + 24)
                            / 6,
                        3,
                        0,
                        &[i, 2 * k + 3 - i],
                        n,
                    )
                }))
                .chain(term(k_ - 2, 0, 1, &[2 * k + 3], n))
                .chain(term(-(2 * k_ - 6), 0, 1, &[1, 2 * k + 2], n))
                .chain(term(2 * k_ - 8, 0, 1, &[2, 2 * k + 1], n))
                .chain(term(-(2 * k_ - 10), 0, 1, &[3, 2 * k], n))
                .chain((4..=k + 1).flat_map(|i| {
                    let i_ = Coefficient::from(i);
                    term(
                        Coefficient::pow(-1, i.into()) * (2 * (k_ - i_) + 3),
                        0,
                        1,
                        &[i, 2 * k + 3 - i],
                        n,
                    )
                }))
                .map(|(m, c)| {
                    (
                        m,
                        (Coefficient::pow(-1, (k + 1).into()) * c) & ((1 << t) - 1),
                    )
                })
                .filter(|(m, c)| *c != 0 && m.zs == m.zs & ((1 << n) - 1))
                .collect(),
        )
    }))
    .collect::<HashMap<_, _>>();

    let mut cache = Default::default();
    let d = u8::try_from(u16::from(n) * u16::from(n + 1) / 2).unwrap();
    for mask in (0..1 << (n - 1)).rev() {
        let cs = (2..=n)
            .filter(|i| mask & (1 << (i - 2)) != 0)
            .collect::<Vec<_>>();
        for d1 in 0..=d - cs.iter().sum::<u8>() {
            if d1 + cs.iter().sum::<u8>() >= d - 3 {
                println!(
                    "z_1^{{{d1}}}c_{{{}}}^*≡{}",
                    cs.iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(","),
                    rational_element(d1, &cs, n, t, &squares, &mut cache)
                        .into_iter()
                        .map(|(m, c)| format!("{c}{m}"))
                        .collect::<Vec<_>>()
                        .join("+"),
                );
            }
        }
    }
}
