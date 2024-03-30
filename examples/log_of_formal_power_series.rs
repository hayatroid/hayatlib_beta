// verification-helper: PROBLEM https://judge.yosupo.jp/problem/log_of_formal_power_series

use ac_library::ModInt998244353 as M;
use hayatlib::polynomial::fps::FPS;
use itertools::Itertools;
use proconio::input;

fn main() {
    input! {
        n: usize,
        a: [M; n],
    }
    let b = FPS::from(a).log(n);
    println!("{}", b.coef.iter().join(" "));
}
