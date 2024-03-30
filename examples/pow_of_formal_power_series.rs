// verification-helper: PROBLEM https://judge.yosupo.jp/problem/pow_of_formal_power_series

use ac_library::ModInt998244353 as M;
use hayatlib::polynomial::fps::FPS;
use itertools::Itertools;
use proconio::input;

fn main() {
    input! {
        n: usize,
        m: usize,
        a: [M; n],
    }
    let b = FPS::from(a).pow(m, n);
    println!("{}", b.coef.iter().join(" "));
}
