// verification-helper: PROBLEM https://judge.yosupo.jp/problem/aplusb

use hayatlib::sample::aplusb;
use proconio::input;

fn main() {
    input! {
        a: u64,
        b: u64,
    }
    println!("{}", aplusb(a, b));
}
