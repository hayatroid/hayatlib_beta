// verification-helper: PROBLEM https://judge.yosupo.jp/problem/composition_of_formal_power_series

use ac_library::ModInt998244353 as M;
use hayatlib::polynomial::FPS;
use itertools::Itertools;
use proconio::input;

fn main() {
    input! {
        n: usize,
        a: [M; n],
        b: [M; n],
    }
    let f = FPS::from(a);
    let g = FPS::from(b);
    println!("{}", f.composition(&g, n).coef.iter().join(" "));
}
