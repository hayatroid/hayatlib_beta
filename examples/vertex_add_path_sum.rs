// verification-helper: PROBLEM https://judge.yosupo.jp/problem/vertex_add_path_sum

use ac_library::{Additive, Segtree};
use hayatlib::data_structure::HLD;
use proconio::input;

fn main() {
    input! {
        n: usize,
        q: usize,
        mut a: [u64; n],
        uv: [(usize, usize); n - 1],
    }
    let mut g = vec![vec![]; n];
    for (u, v) in uv {
        g[u].push(v);
        g[v].push(u);
    }
    let hld = HLD::new(&g, 0);

    let mut seg: Segtree<Additive<u64>> = Segtree::new(n);
    for (i, &a) in a.iter().enumerate() {
        seg.set(hld.pos(i), a);
    }
    
    for _ in 0..q {
        input! {
            op: u8,
        }
        if op == 0 {
            input! {
                p: usize,
                x: u64,
            }
            a[p] += x;
            seg.set(hld.pos(p), a[p]);
        } else {
            input! {
                u: usize,
                v: usize,
            }
            let mut ans = 0;
            let (up, down) = hld.path(u, v);
            for range in up {
                ans += seg.prod(range);
            }
            for range in down {
                ans += seg.prod(range);
            }
            println!("{}", ans);
        }
    }
}
