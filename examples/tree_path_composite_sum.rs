// verification-helper: PROBLEM https://judge.yosupo.jp/problem/tree_path_composite_sum

use ac_library::ModInt998244353 as M;
use hayatlib::tree::{ReRooting, ReRootingMonoid};
use itertools::Itertools;
use proconio::input;

fn main() {
    input! {
        n: usize,
        a: [M; n],
        uvbc: [(usize, usize, M, M); n - 1],
    }
    let mut g = vec![vec![]; n];
    for (u, v, b, c) in uvbc {
        g[u].push((v, (b, c)));
        g[v].push((u, (b, c)));
    }
    let a = a.iter().map(|&a| (a, M::new(1))).collect();
    let dp: ReRooting<X> = ReRooting::new(&g, &a);
    println!("{}", (0..n).map(|i| dp.query(i).0).join(" "));
}

struct X;
impl ReRootingMonoid for X {
    type V = (M, M);
    type E = (M, M);
    fn identity() -> Self::V {
        (M::new(0), M::new(0))
    }
    fn merge(a: Self::V, b: Self::V) -> Self::V {
        (a.0 + b.0, a.1 + b.1)
    }
    fn unmerge(a: Self::V, b: Self::V) -> Self::V {
        (a.0 - b.0, a.1 - b.1)
    }
    fn put_edge(x: Self::V, e: Self::E) -> Self::V {
        (e.0 * x.0 + e.1 * x.1, x.1)
    }
    fn put_vertex(x: Self::V, v: Self::V) -> Self::V {
        (x.0 + v.0, x.1 + v.1)
    }
}
