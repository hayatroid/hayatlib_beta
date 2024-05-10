// verification-helper: PROBLEM https://judge.yosupo.jp/problem/jump_on_tree

use hayatlib::data_structure::HLD;
use proconio::input;

fn main() {
    input! {
        n: usize,
        q: usize,
        ab: [(usize, usize); n - 1],
        sti: [(usize, usize, usize); q],
    }
    let mut g = vec![vec![]; n];
    for (a, b) in ab {
        g[a].push(b);
        g[b].push(a);
    }
    let hld = HLD::new(&g, 0);
    for (s, t, i) in sti {
        if let Some(u) = hld.jump(s, t, i) {
            println!("{}", u);
        } else {
            println!("-1");
        }
    }
}
