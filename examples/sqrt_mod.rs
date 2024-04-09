// verification-helper: PROBLEM https://judge.yosupo.jp/problem/sqrt_mod

use hayatlib::math::sqrt_mod;
use proconio::input;

fn main() {
    input! {
        t: usize,
        yp: [(u32, u32); t],
    }
    for (y, p) in yp {
        if let Some(x) = sqrt_mod(y, p) {
            println!("{}", x);
        } else {
            println!("-1");
        }
    }
}
