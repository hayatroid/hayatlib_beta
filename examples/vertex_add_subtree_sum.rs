// verification-helper: PROBLEM https://judge.yosupo.jp/problem/vertex_add_subtree_sum

use ac_library::{Additive, Segtree};
use hayatlib::data_structure::HLD;
use proconio::input;

fn main() {
    input! {
        n: usize,
        q: usize,
        mut a: [u64; n],
        p: [usize; n - 1],
    }
    let mut g = vec![vec![]; n];
    for (i, &p) in p.iter().enumerate() {
        g[i + 1].push(p);
        g[p].push(i + 1);
    }
    let hld = HLD::new(&g, 0);

    let mut seg: Segtree<Additive<u64>> = Segtree::new(n);
    for (i, &a) in a.iter().enumerate() {
        seg.set(hld.pos[i], a);
    }
    
    for _ in 0..q {
        input! {
            op: u8,
        }
        if op == 0 {
            input! {
                u: usize,
                x: u64,
            }
            a[u] += x;
            seg.set(hld.pos[u], a[u]);
        } else {
            input! {
                u: usize,
            }
            println!("{}", seg.prod(hld.subtree(u)));
        }
    }
}
