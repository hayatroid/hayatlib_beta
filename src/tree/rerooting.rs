//! 全方位木 DP．

pub trait ReRootingMonoid {
    /// 頂点に書かれた値の型．計算結果に加え，計算に必要な補助的な値もここに入れる．
    type V: Copy;
    /// 辺に書かれた値の型．
    type E: Copy;
    fn identity() -> Self::V;
    fn merge(a: Self::V, b: Self::V) -> Self::V;
    fn unmerge(a: Self::V, b: Self::V) -> Self::V;
    fn put_edge(x: Self::V, e: Self::E) -> Self::V;
    fn put_vertex(x: Self::V, v: Self::V) -> Self::V;
}

/// 全方位木 DP．
pub struct ReRooting<M: ReRootingMonoid> {
    g: Vec<Vec<(usize, M::E)>>,
    a: Vec<M::V>,
    root: Vec<M::V>,
    reroot: Vec<M::V>,
}

impl<M: ReRootingMonoid> ReRooting<M> {
    /// 木 $G$，頂点の初期状態の集合 $A$ を受け取り，$\mathrm{Rooting}$，$\mathrm{ReRooting}$ をして返す．
    pub fn new(g: &Vec<Vec<(usize, M::E)>>, a: &Vec<M::V>) -> Self {
        let n = g.len();
        let mut res = Self {
            g: g.clone(),
            a: a.clone(),
            root: vec![M::identity(); n],
            reroot: vec![M::identity(); n],
        };
        res.dfs_root(0, 0);
        res.reroot[0] = res.root[0];
        res.dfs_reroot(0, 0);
        res
    }

    fn dfs_root(&mut self, v: usize, p: usize) {
        let len = self.g[v].len();
        for i in 0..len {
            let (c, e) = self.g[v][i];
            if c == p {
                continue;
            }
            self.dfs_root(c, v);
            // 部分木に辺を加えた「部分木もどき」をマージする
            self.root[v] = M::merge(self.root[v], M::put_edge(self.root[c], e));
        }
        // マージした「部分木もどき」に頂点を付加して部分木にする
        self.root[v] = M::put_vertex(self.root[v], self.a[v]);
    }

    fn dfs_reroot(&mut self, v: usize, p: usize) {
        let len = self.g[v].len();
        for i in 0..len {
            let (c, e) = self.g[v][i];
            if c == p {
                continue;
            }
            // c -> v の辺を v -> c に貼りなおす
            self.root[v] = M::unmerge(self.root[v], M::put_edge(self.root[c], e));
            self.root[c] = M::merge(self.root[c], M::put_edge(self.root[v], e));

            self.reroot[c]  = self.root[c];
            self.dfs_reroot(c, v);

            // 戻す
            self.root[c] = M::unmerge(self.root[c], M::put_edge(self.root[v], e));
            self.root[v] = M::merge(self.root[v], M::put_edge(self.root[c], e));
        }
    }

    /// $G$ 上の点 $v$ を根としたときの計算結果を返す．
    pub fn query(&self, v: usize) -> M::V {
        self.reroot[v]
    }
}
