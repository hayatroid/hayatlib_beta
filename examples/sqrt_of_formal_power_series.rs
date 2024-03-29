// verification-helper: PROBLEM https://judge.yosupo.jp/problem/sqrt_of_formal_power_series

use ac_library::ModInt998244353 as M;
use hayatlib::polynomial::fps::FPS;
use itertools::Itertools;
use proconio::input;

fn main() {
    input! {
        n: usize,
        a: [M; n],
    }
    if let Some(b) = FPS::from(a).sqrt(n) {
        println!("{}", b.coef.iter().join(" "));
    } else {
        println!("-1");
    }
}
