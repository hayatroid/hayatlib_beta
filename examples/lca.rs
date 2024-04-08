// verification-helper: PROBLEM https://judge.yosupo.jp/problem/lca

use hayatlib::data_structure::HLD;
use proconio::input;

fn main() {
    input! {
        n: usize,
        q: usize,
        p: [usize; n - 1],
        uv: [(usize, usize); q],
    }
    let mut g = vec![vec![]; n];
    for (i, &p) in p.iter().enumerate() {
        g[i + 1].push(p);
        g[p].push(i + 1);
    }
    let hld = HLD::new(&g, 0);
    
    for (u, v) in uv {
        println!("{}", hld.lca(u, v));
    }
}
