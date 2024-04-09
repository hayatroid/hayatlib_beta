// verification-helper: PROBLEM https://judge.yosupo.jp/problem/vertex_set_path_composite

use ac_library::{ModInt998244353 as M, Monoid, Segtree};
use hayatlib::data_structure::HLD;
use proconio::input;

fn main() {
    input! {
        n: usize,
        q: usize,
        ab: [(M, M); n],
        uv: [(usize, usize); n - 1],
    }
    let mut g = vec![vec![]; n];
    for (u, v) in uv {
        g[u].push(v);
        g[v].push(u);
    }
    let hld = HLD::new(&g, 0);

    let mut seg: Segtree<X> = Segtree::new(n);
    let mut seg_rev: Segtree<X> = Segtree::new(n);
    for (i, &(a, b)) in ab.iter().enumerate() {
        seg.set(hld.pos(i), (a, b));
        seg_rev.set(n - 1 - hld.pos(i), (a, b));
    }
    
    for _ in 0..q {
        input! {
            op: u8,
        }
        if op == 0 {
            input! {
                p: usize,
                (c, d): (M, M),
            }
            seg.set(hld.pos(p), (c, d));
            seg_rev.set(n - 1 - hld.pos(p), (c, d));
        } else {
            input! {
                u: usize,
                v: usize,
                x: M,
            }
            let mut ans = X::identity();
            let (up, down) = hld.path(u, v);
            for range in up {
                ans = X::binary_operation(&ans, &seg_rev.prod(n - 1 - range.end()..=n - 1 - range.start()));
            }
            for range in down {
                ans = X::binary_operation(&ans, &seg.prod(range));
            }
            println!("{}", ans.0 * x + ans.1);
        }
    }
}

struct X;
impl Monoid for X {
    type S = (M, M);
    fn identity() -> Self::S {
        (M::new(1), M::new(0))
    }
    fn binary_operation(a: &Self::S, b: &Self::S) -> Self::S {
        (
            b.0 * a.0,
            b.0 * a.1 + b.1,
        )
    }
}
